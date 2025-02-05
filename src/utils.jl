
using ConcreteStructs, Distributions

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

nothing