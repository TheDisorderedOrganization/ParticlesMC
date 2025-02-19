struct LAMMPS <: MonteCarlo.Format
    extension::String
    function LAMMPS()
        return new(".lammpstrj")
    end
end

function load_configuration(io, format::LAMMPS; m=-1)
    data = readlines(io)
    timestep_index = findfirst(contains("ITEM: TIMESTEP"), data)
    number_of_atoms_index = findfirst(contains("ITEM: NUMBER OF ATOMS"), data)
    box_bounds_index = findfirst(contains("ITEM: BOX BOUNDS"), data)
    atoms_index = findfirst(contains("ITEM: ATOMS"), data)
    N = parse(Int, data[number_of_atoms_index + 1])
    box_bounds = data[box_bounds_index + 1:box_bounds_index + 3]
    box_bounds = [parse.(Float64, split(elt)) for elt in box_bounds] # Convert row elements to Float64
    box = [elt[2] - elt[1] for elt in box_bounds]
    atom_properties = split(data[atoms_index], " ")[3:end]
    if issubset(["x", "y", "z"], atom_properties)
        d = 3
    elseif issubset(["x", "y"],  atom_properties)
        d = 2
        box = box[1:2]
    else
        error("Unsupported number of dimensions")
    end
    frame = data[atoms_index + 1:atoms_index + N]

    sT = typeof(eval(Meta.parse(split(frame[1], " ")[1])))
    species = Vector{sT}(undef, N)
    position = Vector{SVector{d, Float64}}(undef, N)
    
    for i in eachindex(frame)
        split_line = split(frame[i], " ")
        split_line = filter(x -> x != "", split_line)  
        species[i] = eval(Meta.parse(split_line[1]))
        position[i] = SVector{d}(parse.(Float64, split_line[2:end]))
    end
    return species, position, box, []
end

function MonteCarlo.store_trajectory(io, system::Atoms, t, format::LAMMPS)
    println(io, "ITEM: TIMESTEP")
    println(io, t)
    println(io, "ITEM: NUMBER OF ATOMS")
    println(io, system.N)
    println(io, "ITEM: BOX BOUNDS pp pp pp")
    for i in 1:system.d
        println(io, "0.0 $(system.box[i])")
    end
    if system.d == 2
        println(io, "ITEM: ATOMS type x y")
    elseif system.d == 3
        println(io, "ITEM: ATOMS type x y z")
    end
    for (i, p) in enumerate(system.position)
        println(io, "$(system.species[i])")
        for pos in p
             print(" $pos")
        end
    end
    return nothing
end

function MonteCarlo.store_trajectory(io, system::Molecules, t, format::LAMMPS)
    println(io, "ITEM: TIMESTEP")
    println(io, t)
    println(io, "ITEM: NUMBER OF ATOMS")
    println(io, system.N)
    println(io, "ITEM: BOX BOUNDS pp pp pp")
    for i in 1:system.d
        println(io, "0.0 $(system.box[i])")
    end
    println(io, "ITEM: ATOMS molecule type x y z")
    for (i, p) in enumerate(system.position)
        println(io, "$(system.molecule[i]) $(system.species[i]) $(p[1]) $(p[2]) $(p[3])")
    end
    return nothing
end
