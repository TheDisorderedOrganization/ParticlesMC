struct XYZ <: MonteCarlo.Format
    extension::String
    function XYZ()
        return new(".xyz")
    end
end

function load_configuration(io, format::XYZ; m=-1)
    data = readlines(io)  
    N = parse(Int, data[1])  # Number of atoms or entries
    metadata = split(data[2], " ")  # Metadata split into an array
    
    # Extract cell vector from metadata
    cell_string = replace(metadata[findfirst(startswith("cell:"), metadata)], "cell:" => "")
    cell_vector = parse.(Float64, split(cell_string, ","))
    d = length(cell_vector)
    box = SVector{d}(cell_vector)  # Convert to static vector

    # Determine row index for selected frame
    selrow = m ≥ 0 ? (N + 2) * m - N + 1 : length(data) + m * (N + 2) + 3
    frame = data[selrow:selrow+N-1]

    # Determine species type dynamically
    sT = typeof(eval(Meta.parse(join(split(frame[1], " ")[1:end-d], " "))))
    species = Vector{sT}(undef, N)
    position = Vector{SVector{d, Float64}}(undef, N)

    # Parse species and positions
    for i in eachindex(frame)
        species[i] = eval(Meta.parse(join(split(frame[i], " ")[1:end-d], " ")))
        position[i] = SVector{d}(parse.(Float64, split(frame[i], " ")[end-d+1:end]))
    end
    return species, position, box, metadata
end

function MonteCarlo.store_trajectory(io, system::Atoms, t, format::XYZ)
    println(io, system.N)
    box = replace(replace(string(system.box), r"[\[\]]" => ""), r",\s+" => ",")
    println(io, "step:$t columns:species,position dt:1 cell:$(box) rho:$(system.density) T:$(system.temperature) model:$(system.model.name) potential_energy_per_particle:$(mean(system.local_energy)/2)")
    for (i, p) in enumerate(system.position)
        print(io, "$(system.species[i])")
        for a in 1:system.d
            print(io, " $(p[a])")
        end
        println(io)
    end
    return nothing
end

function MonteCarlo.store_trajectory(io, system::Molecules, t, format::XYZ)
    println(io, system.N)
    box = replace(replace(string(system.box), r"[\[\]]" => ""), r",\s+" => ",")
    println(io, "step:$t columns:molecule,species,position dt:1 cell:$(box) rho:$(system.density) T:$(system.temperature) model:$(system.model.name) potential_energy_per_particle:$(mean(system.local_energy)/2)")
    for (i, p) in enumerate(system.position)
        print(io, "$(system.molecule[i])")
        print(io, "$(system.species[i])")
        for a in 1:system.d
            print(io, " $(p[a])")
        end
        println(io)
    end
    return nothing
end



