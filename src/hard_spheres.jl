mutable struct HardSpheres{D,T<:AbstractFloat,C<:CellList} <: Particles
    position::Vector{SVector{D,T}}
    species::Vector{T}
    density::T
    temperature::T
    d::Int
    N::Int
    box::SVector{D,T}
    cell_list::C
    particle_ids::Base.OneTo{Int}
end

struct HardCore <: Model end

cutoff(si, sj, ::HardCore) = (si + sj) / 2
cutoff2(si, sj, ::HardCore) = (si + sj)^2 / 4

function System(position, species, density::T, temperature::T, ::HardCore; list_type=EmptyList) where {T<:AbstractFloat}
    @assert length(position) == length(species)
    particle_ids = eachindex(position)
    N = length(particle_ids)
    d = length(Array(position)[1])
    box = @SVector fill(T((N / density)^(1 / d)), d)
    maxcut = maximum(species)
    cell_list = list_type(box, maxcut, N)
    system = HardSpheres(position, species, density, temperature, d, N, box, cell_list, particle_ids)
    build_cell_list!(system, system.cell_list)
    return system
end

function Arianna.delta_log_target_density(_::Nothing, overlaps::Bool, ::HardSpheres{D,T,C}) where {D,T,C}
    return overlaps ? typemin(T) : zero(T)
end

function check_overlaps(system, i, ::EmptyList)
    position_i, species_i = system.position[i], system.species[i]
    @inbounds for j in system.particle_ids
        if i != j
            image = nearest_image_distance(position_i, system.position[j], system.box)
            r2 = dot(image, image)
            if r2 < (species_i + system.species[j])^2 / 4
                return true
            end
        end
    end
    return false
end

function check_overlaps(system, i, cell_list::LinkedList)
    # Get cell of particle i
    position_i, species_i = system.position[i], system.species[i]
    mc = get_cell(position_i, cell_list.cell)
    # Scan the neighbourhood of cell mc (including itself)
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(mc2, cell_list.ncells)
        # Scan atoms in cell c2
        j = cell_list.head[c2]
        while (j != -1)
            if i != j
                image = nearest_image_distance(position_i, system.position[j], system.box)
                r2 = dot(image, image)
                if r2 < (species_i + system.species[j])^2 / 4
                    return true
                end
            end
            j = cell_list.list[j]
        end
    end
    return false
end

function check_overlaps(system)
    for i in system.particle_ids
        overlaps = check_overlaps(system, i, EmptyList())
        if overlaps
            return true
        end
    end
    return false
end

function Arianna.perform_action!(system::HardSpheres, action::Displacement)
    system.position[action.i] = system.position[action.i] + action.δ
    c, c2 = old_new_cell(system, action.i, system.cell_list)
    if c != c2
        update_cell_list!(system, action.i, c, c2, system.cell_list)
    end
    overlaps = check_overlaps(system, action.i, system.cell_list)
    return nothing, overlaps
end

function Arianna.revert_action!(system::HardSpheres, action::Displacement)
    system.position[action.i] = system.position[action.i] + action.δ
    c, c2 = old_new_cell(system, action.i, system.cell_list)
    if c != c2
        update_cell_list!(system, action.i, c, c2, system.cell_list)
    end
    return nothing
end

function Arianna.write_system(io, system::HardSpheres)
    println(io, "\tNumber of particles: $(system.N)")
    println(io, "\tDimensions: $(system.d)")
    println(io, "\tCell: $(system.box)")
    println(io, "\tDensity: $(system.density)")
    println(io, "\tTemperature: $(system.temperature)")
    println(io, "\tCell list: " * replace(string(typeof(system.cell_list)), r"\{.*" => ""))
    println(io, "\tModel: HardCore")
    return nothing
end

# NPT ensemble

function Arianna.delta_log_target_density(x₁::Tuple{T,T}, x₂::Tuple{Bool,T}, system::HardSpheres{D,T,C}) where {D,T,C}
    P, V1 = x₁
    overlaps, V2 = x₂
    overlaps && return typemin(T)
    return - P * (V2 - V1) / system.temperature + (system.N + 1) * log(V2 / V1)
end

mutable struct Barostat{T<:AbstractFloat} <: Action
    P::T
    δV::T
end

function new_cell_list(system, rcut, ::EmptyList)
    return EmptyList()
end

function new_cell_list(system, rcut, ::LinkedList)
    return LinkedList(system.box, rcut, system.N)
end

function Arianna.perform_action!(system::Particles, action::Barostat)
    V1 = prod(system.box)
    V2 = V1 + action.δV
    λ = (V2 / V1)^(1 / system.d)
    system.box = system.box .* λ
    system.position .= system.position .* λ
    system.density = system.N / prod(system.box)
    system.cell_list = new_cell_list(system, maximum(system.species), system.cell_list)
    build_cell_list!(system, system.cell_list)
    overlaps = check_overlaps(system)
    return (action.P, V1), (overlaps, V2)
end

function Arianna.invert_action!(action::Barostat, ::Particles)
    action.δV = -action.δV
    return nothing
end

function Arianna.PolicyGuided.reward(action::Barostat, ::Particles)
    return norm(action.δV)^2
end

function Arianna.log_proposal_density(action::Barostat, ::SimpleGaussian, parameters, system::Particles)
    return -action.δV / (2parameters.σ^2) - log(2π * parameters.σ^2) / 2
end

function Arianna.sample_action!(action::Barostat, ::SimpleGaussian, parameters, system::Particles, rng)
    action.δV = rand(rng, Normal(zero(system.density), parameters.σ))
    return nothing
end