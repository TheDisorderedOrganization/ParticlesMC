struct LAMMPS <: MonteCarlo.Format
    extension::String
    function LAMMPS()
        return new(".lammpstrj")
    end
end

function parse_column_string(column_str::AbstractString, ::LAMMPS)
    columns = split(column_str, " ")
    column_info = OrderedDict{String, Vector}() # Use OrderedDict to maintain order
    index = 1
    for column_name in columns
        if column_name == "molecule"
            dimension = 1
            column_info[column_name] = [dimension, index]
        elseif column_name == "type"
            dimension = 1
            column_info["species"] = [dimension, index]
        elseif column_name == "x"
            if issubset(["x", "y", "z"], columns)
                dimension = 3
            elseif issubset(["x", "y"], columns)
                dimension = 2
            column_info["pos"] =  [dimension, index]
            end
        elseif (column_name == "y") ||  (column_name == "y") ||  (column_name == "ITEM:") ||  (column_name == "ATOMS") 
            continue
        else
            error("$column_name is not supported")
        end
        index += 1
    end
    return column_info
end


function get_selrow(::LAMMPS, N, m)
    return m ≥ 0 ? (N + 9) * m - N + 1 : length(data) + m * (N + 9) + 3
end

function get_system_column(::Atoms, ::LAMMPS)
    return ""
end

function get_system_column(::Molecules, ::LAMMPS)
    return "molecule"
end

function read_header(data, format::LAMMPS)
    timestep_index = findfirst(contains("ITEM: TIMESTEP"), data)
    number_of_atoms_index = findfirst(contains("ITEM: NUMBER OF ATOMS"), data)
    box_bounds_index = findfirst(contains("ITEM: BOX BOUNDS"), data)
    columns_index = findfirst(contains("ITEM: ATOMS"), data)
    N = parse(Int, data[number_of_atoms_index + 1])
    box_bounds = data[box_bounds_index + 1:box_bounds_index + 3]
    box_bounds = [parse.(Float64, split(elt)) for elt in box_bounds] # Convert row elements to Float64
    box = [elt[2] - elt[1] for elt in box_bounds]
    column_info = parse_column_string(data[columns_index], format)
    return N, box, column_info, []
end

function write_header(io, system::Particles, t, format::LAMMPS, digits::Integer)
    println(io, "ITEM: TIMESTEP")
    println(io, t)
    println(io, "ITEM: NUMBER OF ATOMS")
    println(io, system.N)
    println(io, "ITEM: BOX BOUNDS pp pp pp")
    for i in 1:system.d
        println(io, "0.0 $(system.box[i])")
    end
    if system.d == 2
        println(io, "-0.1 0.1")
    end
    if system.d == 2
        println(io, "ITEM: ATOMS $(get_system_column(system, format)) type x y")
    elseif system.d == 3
        println(io, "ITEM: ATOMS $(get_system_column(system, format)) type x y z")
    end
    return nothing
end