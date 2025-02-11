
using ConcreteStructs, Distributions, Statistics

function MonteCarlo.delta_log_target_density(e1, e2, system::Particles)
    return -(e2 - e1) ./ system.temperature
end

function nearest_image_distance(xi::T, xj::T, L::T) where {T<:AbstractFloat}
    dx = xi - xj
    return dx - round(dx / L) * L
end

function nearest_image_distance(xi::T, xj::T, box::T) where {T<:AbstractArray}
    return map((xi, xj, L) -> nearest_image_distance(xi, xj, L), xi, xj, box)
end

struct SpeciesList
    sp_ids::Vector{Vector{Int}}
    sp_heads::Vector{Int}
end

function SpeciesList(species)
    available_species = sort(unique(species))
    ids = Vector{Vector{Int}}(undef, length(available_species))
    heads = zeros(Int, length(species))
    for k in eachindex(available_species)
        ids[k] = findall(isequal(available_species[k]), species)
        for i in eachindex(species)
            if species[i] == available_species[k]
                heads[i] = findfirst(isequal(i), ids[k])
            end
        end
    end
    return SpeciesList(ids, heads)
end

function callback_energy(simulation)
    return mean(mean(system.local_energy) for system in simulation.chains) / 2
end

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

function MonteCarlo.store_trajectory(trj, system::Particles, t)
    println(trj, system.N)
    box = replace(replace(string(system.box), r"[\[\]]" => ""), r",\s+" => ",")
    println(trj, "step:$t columns:species,position dt:1 cell:$(box) rho:$(system.density) T:$(system.temperature) model:$(system.model.name) potential_energy_per_particle:$(mean(system.local_energy)/2)")
    for (i, p) in enumerate(system.position)
        print(trj, "$(system.species[i])")
        for a in 1:system.d
            print(trj, " $(p[a])")
        end
        println(trj)
    end
    return nothing
end

function load_configuration(path; m=-1)
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

function volume_sphere(r, d::Int)
    d == 0 && return 1
    d == 1 && return 2r
    return 2π * r^2 * volume_sphere(r, d - 2) / d 
end

function load_init_files(init_path; args=Dict(), tails=".xyz", verbose=false)
    input_files = Vector{String}()
    if endswith(basename(init_path), tails)
        push!(input_files, init_path)
    elseif isdir(init_path)
        for (root, dirs, files) in walkdir(init_path)
            for file in files
                if endswith(file, tails)
                    push!(input_files, joinpath(root, file))
                end
            end
        end
    else
        @error "Unknown format"
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
    initial_temperature_array = [parse(Float64, split(filter(x -> occursin("T:", x), metadata)[1], ":")[2]) for metadata in metadata_array]
    input_models = [split(filter(x -> occursin("model:", x), metadata)[1], ":")[2] for metadata in metadata_array]
    @assert all(isequal(input_models[1]), input_models)
    # Update density, temperature and model if needed
    if haskey(args, "density") && !isnothing(args["density"])
        λs = (initial_density_array ./ args["density"]) .^ (1 / d)
        initial_density_array .= args["density"]
        initial_position_array .= [X .* λ for (X, λ) in zip(initial_position_array, λs)]
        initial_box_array .= [box .* λ for (box, λ) in zip(initial_box_array, λs)]
    end
    if haskey(args, "temperature") && !isnothing(args["temperature"])
        initial_temperature_array .= args["temperature"]
    end
    if haskey(args, "model") && !isnothing(args["model"])
        input_models .= args["model"]
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
    # Create independent chains
    chains = [System(initial_position_array[k], initial_species_array[k], initial_density_array[k], initial_temperature_array[k], model, list_type=list_type) for k in eachindex(initial_position_array)]
    verbose && println("$(length(chains)) chains created")
    return chains
end

nothing