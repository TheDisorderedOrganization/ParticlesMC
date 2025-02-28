module IO

using ..ParticlesMC: Particles, Atoms, Molecules, System
using ..ParticlesMC: fold_back, cutoff, volume_sphere
using ..ParticlesMC: EmptyList, LinkedList
using ..ParticlesMC: Model, GeneralKG, JBB, BHHP, SoftSpheres, KobAndersen
using Arianna
using Distributions, LinearAlgebra, StaticArrays, Printf
using DataStructures: OrderedDict
export XYZ, EXYZ, LAMMPS
export load_configuration, load_chains

include("xyz.jl")
include("exyz.jl")
include("lammps.jl")

function Arianna.write_system(io, system::Particles)
    println(io, "\tNumber of particles: $(system.N)")
    println(io, "\tDimensions: $(system.d)")
    println(io, "\tCell: $(system.box)")
    println(io, "\tDensity: $(system.density)")
    println(io, "\tTemperature: $(system.temperature)")
    println(io, "\tCell list: " * replace(string(typeof(system.cell_list)), r"\{.*" => ""))
    println(io, "\tModel: $(system.model.name)")
    return nothing
end

function load_configuration(filename::String)
    io = open(filename, "r")  # Open file as IOStream
    if endswith(filename, ".xyz")
        return load_configuration(io, XYZ())
    elseif endswith(filename, ".exyz")
        return load_configuration(io, EXYZ())
    elseif endswith(filename, ".lmp") || endswith(filename, ".lammpstrj") || endswith(filename, ".lammps")
        return load_configuration(io, LAMMPS())
    else
        error("Unsupported file format: $filename")
    end
end

function load_configuration(io, format::Arianna.Format; m=1)
    data = readlines(io)
    N, box, column_info, metadata = read_header(data, format)
    selrow = get_selrow(format, N, m)
    frame = data[selrow:selrow+N-1]
    bool_molecule = "molecule" in keys(column_info)
    bool_species = "species" in keys(column_info)
    if bool_molecule
        if m != 1
            error("For molecular systems the frame index has to be equal to 1")
        end
        molecule_d, molecule_index = column_info["molecule"]
        if molecule_d != 1
            error("molecule dimension must be 1")
        end
        molecule = Vector{Int}(undef, N)
        btype, bond = read_bonds(data, N, format)
    end
    if bool_species
        species_d, species_index = column_info["species"]
        if species_d != 1
            error("Species dimension must be 1")
        end
        sT = typeof(eval(Meta.parse(split(frame[1], " ")[1])))
        species = Vector{sT}(undef, N)
    else
        species = ones(Int, N)
    end
    if "pos" in keys(column_info)
        pos_d, pos_index = column_info["pos"]
        position = Vector{SVector{pos_d, Float64}}(undef, N)
    else
        missing_key_error("pos")
    end
    if pos_d < length(box)
        box = box[1:pos_d]
    end
    for i in eachindex(frame)
        split_line = split(frame[i], " ")
        if bool_species
            species[i] = eval(Meta.parse.(split_line[species_index]))
        end
        if bool_molecule
            molecule[i] = parse.(Int64, split_line[molecule_index])
        end
        position[i] = SVector{pos_d}(parse.(Float64, split_line[pos_index:pos_index+pos_d-1]))
    end
    config_dict = Dict( :N => N,
                        :d => pos_d,
                        :box => box,
                        :species => species,
                        :position => position,
                        :metadata => metadata
    )
    if bool_molecule
        config_dict[:molecule] = molecule
        config_dict[:btype] = btype
        config_dict[:bond] = bond
    end
    return config_dict
end

function read_bonds(data, N, format::Arianna.Format)
    selrow = get_selrow(format, N, 1)
    bonds_data = data[N+selrow:end]
    
    if length(bonds_data) == 0
        error("No bonds found in the file")
    else
        N_bonds, column_info = read_bonds_header(bonds_data, format)
    end
    bool_btype = "btype" in keys(column_info)
    bool_bond = "bond" in keys(column_info)
    if bool_btype
        btype_d, btype_index = column_info["btype"]
        if btype_d != 1
            error("Bond type dimension must be 1")
        end
        btype = fill(Vector{Int}(), N)
    end
    if bool_bond
        bond_d, bond_index = column_info["bond"]
        if bond_d != 2
            error("Bond dimension must be 2. Found $bond_d.")
        end
        row_bonds = get_row_bonds(selrow, N, format)
        bond = [Vector{Int}() for _ in 1:N]
        for i in 1:N_bonds
            atom_i, atom_j = parse.(Int, split(bonds_data[row_bonds + i], " ")[bond_index:bond_index+1])
            push!(bond[atom_i], atom_j)
            push!(bond[atom_j], atom_i)
            if bool_btype
                btype_ij = parse.(Int, split(bonds_data[row_bonds + i], " ")[btype_index])
            else
                btype_ij = 1
            end
            push!(btype[atom_i], btype_ij)
            push!(btype[atom_j], btype_ij)
        end
    else
        error("Bond array is not written in the $format file")
    end
    return btype, bond
end

function missing_key_error(key)
    error(error("$key array has not been found in metadata or is not defined. Define the $key in the args Dict"))
end

function broadcast_dict(dicts, key)
    return [dict[key] for dict in dicts]
end

function load_chains(init_path; args=Dict(), verbose=false)
    input_files = Vector{String}()
    if isfile(init_path)
        push!(input_files, init_path)
    elseif isdir(init_path)
        for (root, dirs, files) in walkdir(init_path)
            for file in files
                push!(input_files, joinpath(root, file))
            end
        end
    end
    verbose && println("Processing $(length(input_files)) configuration file(s)")
    verbose && @show input_files
    config_dict = load_configuration.(input_files)
    initial_species_array = broadcast_dict(config_dict, :species)
    initial_position_array =  broadcast_dict(config_dict, :position)
    initial_box_array =  broadcast_dict(config_dict, :box)
    metadata_array = broadcast_dict(config_dict, :metadata)
    N, d = config_dict[1][:N], config_dict[1][:d]
    @assert all(isequal(N), length.(initial_position_array))
    @assert all(isequal(d), vcat([length.(X) for X in initial_position_array]...))
    initial_density_array = length.(initial_position_array) ./ prod.(initial_box_array)
    if length(metadata_array) > 1
        initial_temperature_array = [parse(Float64, split(filter(x -> occursin("T:", x), metadata)[1], ":")[2]) for metadata in metadata_array]
        input_models = [split(filter(x -> occursin("model:", x), metadata)[1], ":")[2] for metadata in metadata_array]
        @assert all(isequal(input_models[1]), input_models)
    else
        initial_temperature_array = nothing
        input_models = nothing
    end
    # Update density, temperature and model if needed
    if haskey(args, "density") && !isnothing(args["density"])
        λs = (initial_density_array ./ args["density"]) .^ (1 / d)
        initial_density_array .= args["density"]
        initial_position_array .= [X .* λ for (X, λ) in zip(initial_position_array, λs)]
        initial_box_array .= [box .* λ for (box, λ) in zip(initial_box_array, λs)]
    end
    if haskey(args, "temperature") && !isnothing(args["temperature"])
        initial_temperature_array = args["temperature"]
    elseif isnothing(initial_temperature_array)
        missing_key_error("temperature")
    end
    if haskey(args, "model") && !isnothing(args["model"])
        input_models = args["model"]

    elseif isnothing(input_models)
        missing_key_error("model")
    end
    # Fold back into the box
    initial_position_array .= [[fold_back(x, box) for x in X] for (X, box) in zip(initial_position_array, initial_box_array)]
    # Parse model
    if occursin(r"\(", input_models[1]) && occursin(r"\)", input_models[1])
        model = eval(Meta.parse(input_models[1]))  # Parse the string if it has parentheses
    else
        model =  eval(Meta.parse(input_models[1] * "()"))  # Else, append () and evaluate
    end
    @assert isa(model, Model)
    # Copy configurations nsim times (replicas)
    if haskey(args, "nsim") && !isnothing(args["nsim"]) && args["nsim"] > 1
        nsim = args["nsim"]
        verbose && println("Generating $nsim replicas per input file")
        initial_position_array = vcat([[copy(x) for _ in 1:nsim] for x in initial_position_array]...)
        initial_species_array = vcat([[copy(x) for _ in 1:nsim] for x in initial_species_array]...)
        initial_density_array = vcat([[copy(x) for _ in 1:nsim] for x in initial_density_array]...)
        initial_temperature_array = vcat([[copy(x) for _ in 1:nsim] for x in initial_temperature_array]...)
    end
    # Handle cell list (this is classy)
    available_species = unique(vcat(initial_species_array...))
    maxcut = maximum([cutoff(spi, spj, model) for spi in available_species for spj in available_species])
    Z = mean(initial_density_array) * volume_sphere(maxcut, d)
    list_type = Z / N < 0.1 ? LinkedList : EmptyList
    if haskey(args, "list_type") && !isnothing(args["list_type"])
        list_type = eval(Meta.parse(args["list_type"]))
    end
    verbose && println("Using $list_type as cell list type")
    # Create independent chains
    bool_molecule = :molecule in keys(config_dict[1])
    if bool_molecule
        initial_molecule_array = broadcast_dict(config_dict, :molecule)
        initial_bond_array = broadcast_dict(config_dict, :bond)
        initial_btype_array = broadcast_dict(config_dict, :btype)
        chains = [System(initial_position_array[k], initial_species_array[k], initial_molecule_array[k], initial_density_array[k], initial_temperature_array[k], model, initial_bond_array[k], list_type=list_type) for k in eachindex(initial_position_array)]
    else
        chains = [System(initial_position_array[k], initial_species_array[k], initial_density_array[k], initial_temperature_array[k], model, list_type=list_type) for k in eachindex(initial_position_array)]
    end
    verbose && println("$(length(chains)) chains created")
    return chains
end

function formatted_string(num::Real, digits::Integer)
    fmtstr = "%." * string(digits) * "f"
    fmt = Printf.Format(fmtstr)
    return Printf.format(fmt, num)
end

function write_position(io, position, digits::Int)
    for position_i in position 
        formatted_position_i = formatted_string(position_i, digits)
        print(io, " ")
        print(io, formatted_position_i)
    end
    println(io)
    return nothing
end

function Arianna.store_trajectory(io, system::Atoms, t, format::Arianna.Format; digits::Integer=6)
    write_header(io, system, t, format, digits)
    for (species, position) in zip(system.species, system.position)
        print(io, "$species")
        write_position(io, position, digits)
    end
    return nothing
end

function Arianna.store_trajectory(io, system::Molecules, t, format::Arianna.Format; digits::Integer=6)
    write_header(io, system, t, format, digits)
    for (molecule, species, position) in zip(system.molecule, system.species, system.position)
        print(io, "$molecule $species")
        write_position(io, position, digits)
    end
    return nothing
end


end # module