struct XYZ <: Arianna.Format
    extension::String
    function XYZ()
        return new(".xyz")
    end
end

function get_selrow(::XYZ, N, m)
    return m ≥ 0 ? (N + 2) * m - N + 1 : length(data) + m * (N + 2) + 3
end

function parse_column_string(column_str::AbstractString, ::XYZ; d::Int=3)
    columns = split(column_str, ",")
    column_info = OrderedDict{String, Vector}() # Use OrderedDict to maintain order
    index = 1
    for column_name in columns
        if column_name == "molecule"
            dimension = 1
            column_info[column_name] = [dimension, index]
        elseif column_name == "species"
            dimension = 1
            column_info[column_name] = [dimension, index]
        elseif column_name == "position"
            column_info["pos"] =  [d, index]
        elseif column_name == "bond"
            dimension = 2
            column_info["bond"] =  [2, index]
        elseif column_name == "btype"
            dimension = 1
            column_info[column_name] = [dimension, index]
        else
            error("$column_name is not supported")
        end
        index += 1
    end
    return column_info
end

function read_header(data, format::XYZ)
    N = parse(Int, data[1])  # Number of atoms or entries
    metadata = split(data[2], " ")  # Metadata split into an array
    
    # Extract cell vector from metadata
    cell_str = replace(metadata[findfirst(startswith("cell:"), metadata)], "cell:" => "")
    cell_vector = parse.(Float64, split(cell_str, ","))
    d = length(cell_vector)
    box = SVector{d}(cell_vector)
    column_str = replace(metadata[findfirst(startswith("columns:"), metadata)], "columns:" => "")
    column_info = parse_column_string(column_str, format; d=d)
    return N, box, column_info, []
end

function get_system_column(::Atoms, ::XYZ)
    return ""
end

function get_system_column(::Molecules, ::XYZ)
    return "molecule,"
end

function get_row_bonds(selrow, N, ::XYZ)
    return 2
end

function read_bonds_header(bonds, format::XYZ)
    N_bonds = parse(Int, bonds[1])
    column_line = split(bonds[2], " ")
    column_line isa String ? [column_line] : column_line
    column_str = replace(column_line[findfirst(startswith("columns:"), column_line)], "columns:" => "")
    column_info = parse_column_string(column_str, format)
    return N_bonds, column_info
end


function write_header(io, system::Particles, t, format::XYZ, digits::Integer)
    println(io, length(system))
    box = replace(replace(string(system.box), r"[\[\]]" => ""), r",\s+" => ",")
    println(io, "step:$t columns:$(get_system_column(system, format))species,position dt:1 cell:$(box) rho:$(system.density) T:$(system.temperature) potential_energy_per_particle:$(mean(system.local_energy)/2)")
    return nothing
end

