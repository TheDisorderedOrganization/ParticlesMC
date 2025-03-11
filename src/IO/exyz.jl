struct EXYZ <: Arianna.Format
    extension::String
    function EXYZ()
        return new(".exyz")
    end
end

function parse_column_string(column_str::AbstractString, ::EXYZ)
    columns = split(column_str, ":")
    column_info = OrderedDict{String, Vector}() # Use OrderedDict to maintain order
    i, index = 1, 1
    types = ["S", "I", "R"]
    while i <= length(columns)
        if i + 2 <= length(columns) && (columns[i+1] ∈ types)
          column_name = columns[i]
          dimension = parse(Int, columns[i + 2])
          column_info[column_name] = [dimension, index]
          index += dimension
          i += 3 # Skip data type and dimension
        else
          i += 1
        end
    end

    return column_info
end

function read_header(data, format::EXYZ)
    N = parse(Int, data[1])
    metadata_line = data[2]
    mat = match(r"Lattice=\"(.*?)\"", metadata_line)
    if mat === nothing
        error("Invalid Lattice line format")
    end
    lattice_str = mat.captures[1]
    # Convert lattice string to vector of floats
    lattice_values = map(x -> parse(Float64, x), split(lattice_str))
    # Construct lattice matrix
    if length(lattice_values) != 9
        error("Lattice matrix must have 9 elements")
    end
    lattice_matrix = reshape(lattice_values, 3, 3)
    box = lattice_matrix[diagind(lattice_matrix)]
    column_match = match(r"Properties=(.*)", metadata_line)
    column_str = mat === nothing ? nothing : column_match.captures[1]
    column_info =  parse_column_string(column_str, format)
    return N, box, column_info, []
end

function get_selrow(::EXYZ, N, m)
    return m ≥ 0 ? (N + 2) * m - N + 1 : length(data) + m * (N + 2) + 3
end

function compute_box_str(box, ::EXYZ)
    if length(box) == 2
        return "$(box[1]) 0.0 0.0 0.0 $(box[2]) 0.0 0.0 0.0 0.0"
    elseif length(box) == 3
        return "$(box[1]) 0.0 0.0 0.0 $(box[2]) 0.0 0.0 0.0 $(box[3])"
    else
        throw(ArgumentError("Box vector must have 2 or 3 elements."))
    end
end


function get_system_column(::Atoms, ::EXYZ)
    return ""
end

function get_system_column(::Molecules, ::EXYZ)
    return "molecule:I:1"
end

function get_row_bonds(selrow, N, ::EXYZ)
    return 2
end

function read_bonds_header(bonds, format::EXYZ)
    N_bonds = parse(Int, bonds[1])
    metadata_line = bonds[2]
    column_match = match(r"Properties=(.*)", metadata_line)
    column_str = column_match === nothing ? nothing : column_match.captures[1]
    column_info =  parse_column_string(column_str, format)
    return N_bonds, column_info
end

function write_header(io, system::Particles, t, format::EXYZ, digits::Integer)
    println(io, length(system))
    box_str = compute_box_str(system.box, format)
    println(io, "Lattice=\"$box_str\" Properties=$(get_system_column(system, format)):species:S:1:pos:R:$(system.d) Time=$t")
    return nothing
end