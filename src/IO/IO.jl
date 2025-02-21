module IO

using ..ParticlesMC: Particles, Atoms, Molecules, System
using ..ParticlesMC: fold_back, cutoff, volume_sphere
using ..ParticlesMC: EmptyList, LinkedList
using ..ParticlesMC: Model, GeneralKG, JBB, BHHP, SoftSpheres, KobAndersen
using MonteCarlo
using Distributions, LinearAlgebra, StaticArrays, Printf
using DataStructures: OrderedDict
export XYZ, EXYZ, LAMMPS
export load_configuration, load_chains

include("xyz.jl")
include("exyz.jl")
include("lammps.jl")

function MonteCarlo.write_system(io, system::Particles)
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

function missing_key_error(key)
    error(error("$key array has not been found in metadata or is not defined. Define the $key in the args Dict"))
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
    input_data = load_configuration.(input_files)
    initial_species_array = getindex.(input_data, 1)
    initial_position_array = getindex.(input_data, 2)
    initial_box_array = getindex.(input_data, 3)
    metadata_array = getindex.(input_data, 4)
    N = first(length.(initial_position_array))
    d = first(length.(first(initial_position_array)))
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
    model = eval(Meta.parse(input_models[1] * "()"))
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
    chains = [System(initial_position_array[k], initial_species_array[k], initial_density_array[k], initial_temperature_array[k], model, list_type=list_type) for k in eachindex(initial_position_array)]
    verbose && println("$(length(chains)) chains created")
    return chains
end

function formatted_string(num::Real, digits::Integer)
    fmtstr = "%." * string(digits) * "f"
    fmt = Printf.Format(fmtstr)
    return Printf.format(fmt, num)
end

function MonteCarlo.store_trajectory(io, system::Atoms, t, format::MonteCarlo.Format, digits::Integer=6)
    write_header(io, system, t, format, digits)
    for (species, position) in zip(system.species, system.position)
        print(io, "$species")
        for position_i in position 
            formatted_position_i = formatted_string(position_i, digits)
            print(io, " ")
            print(io, formatted_position_i)
        end
        println(io)
    end
    return nothing
end

function MonteCarlo.store_trajectory(io, system::Molecules, t, format::MonteCarlo.Format, digits::Integer=6)
    write_header(io, system, t, format, digits)
    for (molecule, species, position) in zip(system.molecule, system.species, system.position)
        print(io, "$molecule $species")
        for position_i in position 
            formatted_position_i = formatted_string(position_i, digits)
            print(io, " ")
            print(io, formatted_position_i)
        end
        println(io)
    end
    return nothing
end


end