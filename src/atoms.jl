struct Atoms{D, V<:AbstractVector, C<:NeighbourList, H, T<:AbstractFloat, SM<:AbstractArray} <: Particles
    position::Vector{SVector{D,T}}
    species::V
    density::T
    energy::Vector{T}
    temperature::T
    model_matrix::SM
    N::Int
    d::Int
    box::SVector{D,T}
    neighbour_list::C
    species_list::H
end


function System(position, species, density::T, temperature::T, model_matrix; list_type=EmptyList) where{T<:AbstractFloat}
    @assert length(position) == length(species)
    N = length(position)
    d = length(Array(position)[1])
    energy = Vector{T}(undef, 1)
    box = @SVector fill(T((N / density)^(1 / d)), d)
    maxcut = maximum([model.rcut for model in model_matrix])
    neighbour_list = list_type(box, maxcut, N)
    species_list = isa(species[1], Integer) ? SpeciesList(species) : nothing
    system = Atoms(position, species, density, energy, temperature, model_matrix, N, d, box, neighbour_list, species_list)
    build_neighbour_list!(system)
    local_energy = [compute_energy_particle(system, i, neighbour_list) for i in eachindex(position)]
    energy = sum(local_energy) / 2
    if isinf(energy) || isnan(energy)
        error("Initial configuration has infinite or NaN energy.")
    end
    system.energy[1] = energy
    return system
end

function compute_energy_ij(system::Atoms, i, j, position_i)
    # Early return using && for short-circuit evaluation
    i == j && return zero(typeof(system.density))
    # If i != j, compute energy directly
    model_ij = get_model(system, i, j)
    position_j = get_position(system, j)
    r2 = nearest_image_distance_squared(position_i, position_j, get_box(system))
    r2 > cutoff2(model_ij) && return zero(typeof(system.density))
    return potential(r2, model_ij)
end


function compute_energy_particle(system::Atoms, i)
    return compute_energy_particle(system, i, system.neighbour_list)
end

function compute_energy_particle(system::Atoms, i, ::EmptyList)
    energy_i = zero(typeof(system.density))
    position_i = get_position(system, i)
    for (j, _) in enumerate(system)
        energy_i += compute_energy_ij(system, i, j, position_i)
    end
    return energy_i
end


function compute_energy_particle(system::Atoms, i, neighbour_list::CellList)
    energy_i = zero(typeof(system.density))
    position_i = get_position(system, i)
    c = get_cell_index(position_i, neighbour_list)
    neighbour_cells = neighbour_list.neighbour_cells[c]

    # Scan the neighbourhood of cell mc (including itself)
    @inbounds for c2 in neighbour_cells
        # Scan atoms in cell c2
        neighbours = neighbour_list.cells[c2]
        @inbounds for j in neighbours
            energy_i += compute_energy_ij(system, i, j, position_i)
        end
    end
    return energy_i
end

function compute_energy_particle(system::Atoms, i, neighbour_list::LinkedList)
    energy_i = zero(typeof(system.density))
    # Get cell of particle i
    position_i = get_position(system, i)
    c = get_cell_index(position_i, neighbour_list)
    neighbour_cells = neighbour_list.neighbour_cells[c]
    # Scan the neighbourhood of cell mc (including itself)
    @inbounds for c2 in neighbour_cells
        # Scan atoms in cell c2
        j = neighbour_list.head[c2]
        while (j != -1)
            energy_ij = compute_energy_ij(system, i, j, position_i)
            energy_i += energy_ij
            j = neighbour_list.list[j]
        end
    end
    return energy_i
end