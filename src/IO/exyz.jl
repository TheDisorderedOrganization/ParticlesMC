struct EXYZ <: MonteCarlo.Format
    extension::String
    function EXYZ()
        return new(".exyz")
    end
end

#=
function load_configuration(path, format::EXYZ; m=-1)
    data = readlines(path)
    N = parse(Int, data[1])
    metadata = split(data[2], " ")
    cell_string = replace(metadata[findfirst(startswith("cell:"), metadata)], "cell:" => "")
    cell_vector = parse.(Float64, split(cell_string, ","))
    d = length(cell_vector)
    box = SVector{d}(cell_vector)
    selrow = m ≥ 0 ? (N + 2) * m - N + 1 : length(data) + m * (N + 2) + 3
    frame = data[selrow:selrow+N-1]
    sT = typeof(eval(Meta.parse(join(split(frame[1], " ")[1:end-d], " "))))
    species = Vector{sT}(undef, N)
    position = Vector{SVector{d}{Float64}}(undef, N)
    for i in eachindex(frame)
        species[i] = eval(Meta.parse(join(split(frame[i], " ")[1:end-d], " ")))
        position[i] = SVector{d}(parse.(Float64, split(frame[i], " ")[end-d+1:end]))
    end
    return species, position, box, metadata
end

function MonteCarlo.store_trajectory(io, system::Particles, t, format::EXYZ)
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
=#