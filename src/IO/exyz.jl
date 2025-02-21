struct EXYZ <: MonteCarlo.Format
    extension::String
    function EXYZ()
        return new(".exyz")
    end
end

function parse_properties_string(properties_str::AbstractString)
    properties = split(properties_str, ":")
    column_info = OrderedDict{String, Vector}() # Use OrderedDict to maintain order
    i, index = 1, 1
    types = ["S", "I", "R"]
    while i <= length(properties)
        if i + 2 <= length(properties) && (properties[i+1] ∈ types)
          column_name = properties[i]
          dimension = parse(Int, properties[i + 2])
          column_info[column_name] = [dimension, index]
          index += dimension
          i += 3 # Skip data type and dimension
        else
          i += 1
        end
    end

    return column_info
end

function load_configuration(io, format::EXYZ; m=1)
    data = readlines(io)
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

    selrow = m ≥ 0 ? (N + 2) * m - N + 1 : length(data) + m * (N + 2) + 3
    frame = data[selrow:selrow+N-1]

    sT = typeof(eval(Meta.parse(split(frame[1], " ")[1])))
    mat = match(r"Properties=(.*)", metadata_line)
    properties_str = mat === nothing ? nothing : mat.captures[1]
    properties =  parse_properties_string(properties_str)
    if "species" in keys(properties)
        species_d, species_index = properties["species"]
        if species_d > 1
            error("Species dimension must be 1")
        end
        species = Vector{sT}(undef, N)
    end
    if "pos" in keys(properties)
        pos_d, pos_index = properties["pos"]
        position = Vector{SVector{pos_d, Float64}}(undef, N)
    else
        missing_key_error("pos")
    end
    if pos_d < length(box)
        box = box[1:pos_d]
    end
    for i in eachindex(frame)
        split_line = split(frame[i], " ")
        species[i] = eval(Meta.parse.(split_line[species_index:species_index+species_d-1])[1])
        position[i] = SVector{pos_d}(parse.(Float64, split_line[pos_index:pos_index+pos_d-1]))
    end
    return species, position, box, []
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


function get_system_column(system::Atoms, ::EXYZ)
    return ""
end

function get_system_column(system::Molecules, ::EXYZ)
    return "molecule:I:1"
end

function write_header(io, system::Particles, t, format::EXYZ, digits::Integer)
    println(io, system.N)
    box_str = compute_box_str(system.box, format)
    println(io, "Lattice=\"$box_str\" Properties=:$(get_system_column(system, format)):species:I:1:pos:R:$(system.d) Time=$t")
    return nothing
end