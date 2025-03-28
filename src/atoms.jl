struct Atoms{D, Q, V<:AbstractVector, C<:NeighbourList, H, T<:AbstractFloat, SM<:AbstractArray} <: Particles
    position::Vector{SVector{D,T}}
    species::V
    density::T
    energy::Vector{T}
    temperature::T
    model_matrix::SM
    N::Int
    d::Int
    box::SVector{D,T}
    local_energy::MVector{Q, T}
    cell_list::C
    species_list::H
end


function System(position, species, density::T, temperature::T, model_matrix; list_type=EmptyList) where{T<:AbstractFloat}
    @assert length(position) == length(species)
    N = length(position)
    d = length(Array(position)[1])
    energy = Vector{T}(undef, 1)
    box = @SVector fill(T((N / density)^(1 / d)), d)
    local_energy = MVector{N, T}(zeros(T, N))
    indices = MVector{N, Bool}(falses(N))
    maxcut = maximum([model.rcut for model in model_matrix])
    cell_list = list_type(box, maxcut, N)
    species_list = isa(species[1], Integer) ? SpeciesList(species) : nothing
    system = Atoms(position, species, density, energy, temperature, model_matrix, N, d, box, local_energy, cell_list, species_list)
    build_cell_list!(system)
    system.local_energy .= [compute_local_energy(system, i) for i in eachindex(system)]
    system.energy[1] = sum(system.local_energy) / 2
    return system
end

function check_compute_energy_ij(system::Atoms, i, j, position_i)
    # Early return using && for short-circuit evaluation
    i == j && return zero(typeof(system.density))
    # If i != j, compute energy directly
    model_ij = get_model(system, i, j)
    position_j = get_position(system, j)
    return compute_energy_ij(system, position_i, position_j, model_ij)
end

function compute_energy_ij(system::Atoms, position_i, position_j, model_ij::Model)
    r2 = nearest_image_distance_squared(position_i, position_j, get_box(system))
    r2 > cutoff2(model_ij) && return zero(typeof(system.density))
    return potential(r2, model_ij)
end

function compute_local_energy(system::Atoms, i)
    return compute_local_energy(system, i, system.cell_list)
end

function compute_local_energy(system::Atoms, i, ::EmptyList)
    energy_i = zero(typeof(system.density))
    position_i = get_position(system, i)
    for (j, _) in enumerate(system)
        energy_i += check_compute_energy_ij(system, i, j, position_i)
    end
    return energy_i
end


function destroy_particle!(system::Atoms, i, ::EmptyList)
    energy_i = zero(typeof(system.density))
    position_i = get_position(system, i)
    # Loop over particles
    @inbounds for (j, _) in enumerate(system)
        energy_ij = check_compute_energy_ij(system, i, j, position_i)
        energy_i += energy_ij
    end
    return energy_i
end

function create_particle!(system::Atoms, i, ::EmptyList)
    energy_i = zero(typeof(system.density))
    position_i = get_position(system, i)
    # Loop over particles
    @inbounds for (j, position_j) in enumerate(system)
        energy_ij = check_compute_energy_ij(system, i, j, position_i)
        energy_i += energy_ij
    end
    return energy_i
end

function compute_local_energy(system::Atoms, i, cell_list::CellList)
    energy_i = zero(typeof(system.density))
    position_i = get_position(system, i)
    mc = get_cell(position_i, cell_list)
    # Scan the neighbourhood of cell mc (including itself)
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(cell_list, mc2)
        # Scan atoms in cell c2
        @inbounds for j in cell_list.cells[c2]
            energy_i += check_compute_energy_ij(system, i, j, position_i)
        end
    end
    return energy_i
end

# With linked list
function destroy_particle!(system::Atoms, i, cell_list::CellList)
    # Get cell of particle i
    energy_i = zero(typeof(system.density))
    position_i = get_position(system, i)
    mc = get_cell(position_i, cell_list)
    # Scan the neighbourhood of cell mc (including itself)
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(cell_list, mc2)
        # Scan atoms in cell c2
        neighbours = cell_list.cells[c2]
        @inbounds for j in neighbours
            energy_ij = check_compute_energy_ij(system, i, j, position_i)
            energy_i += energy_ij
        end
    end
    return energy_i
end

function create_particle!(system::Atoms, i, cell_list::CellList)
    energy_i = zero(typeof(system.density))
    # Get cell of particle i
    position_i = get_position(system, i)
    mc = get_cell(position_i, cell_list)
    # Scan the neighbourhood of cell mc (including itself)
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(cell_list, mc2)
        # Scan atoms in cell c2
        neighbours = cell_list.cells[c2]
        @inbounds for j in neighbours
            energy_ij = check_compute_energy_ij(system, i, j, position_i)
            energy_i += energy_ij
        end
    end
    return energy_i
end


function compute_local_energy(system::Atoms, i, cell_list::LinkedList)
    energy_i = zero(typeof(system.density))
    # Get cell of particle i
    position_i = get_position(system, i)
    mc = get_cell(position_i, cell_list)
    # Scan the neighbourhood of cell mc (including itself)
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(cell_list, mc2)
        # Scan atoms in cell c2
        j = cell_list.head[c2]
        while (j != -1)
            energy_ij = check_compute_energy_ij(system, i, j, position_i)
            energy_i += energy_ij
            j = cell_list.list[j]
        end
    end
    return energy_i
end


function destroy_particle!(system::Atoms, i, cell_list::LinkedList)
    energy_i = zero(typeof(system.density))
    # Get cell of particle i
    position_i = get_position(system, i)
    mc = get_cell(position_i, cell_list)
    # Scan the neighbourhood of cell mc (including itself)
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(cell_list, mc2)
        # Scan atoms in cell c2
        j = cell_list.head[c2]
        while (j != -1)
            energy_ij = check_compute_energy_ij(system, i, j, position_i)
            energy_i += energy_ij
            j = cell_list.list[j]
        end
    end
    return energy_i
end

function create_particle!(system::Atoms, i, cell_list::LinkedList)
    energy_i = zero(typeof(system.density))
    # Get cell of particle i
    position_i = get_position(system, i)
    mc = get_cell(position_i, cell_list)
    # Scan the neighbourhood of cell mc (including itself)
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(cell_list, mc2)
        # Scan atoms in cell c2
        j = cell_list.head[c2]
        while (j != -1)
            energy_ij = check_compute_energy_ij(system, i, j, position_i)
            energy_i += energy_ij
            j = cell_list.list[j]
        end
    end
    return energy_i
end