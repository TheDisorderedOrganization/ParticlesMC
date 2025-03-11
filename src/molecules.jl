struct Molecules{D,  VS<:AbstractVector, M<:Model, C<:NeighbourList, T<:AbstractFloat} <: Particles
    position::Vector{SVector{D,T}}
    species::VS
    molecule::VS
    molecule_species::VS
    density::T
    temperature::T
    energy::Vector{T}
    model_matrix::Matrix{M}
    d::Int
    N::Int
    box::SVector{D,T}
    local_energy::Vector{T}
    cell_list::C
    bonds::Vector{Vector{Int}}
end


function System(position, species, molecule, density::T, temperature::T, model_matrix::Matrix{M}, bonds; molecule_species=nothing, list_type=EmptyList) where {T<:AbstractFloat, M<:Model}
    @assert length(position) == length(species)
    particle_ids = eachindex(position)
    N = length(particle_ids)
    molecule_species = something(molecule_species, ones(Int, N))
    d = length(Array(position)[1])
    box = @SVector fill(T((N / density)^(1 / d)), d)
    e_locals = zeros(T, N)
    energy = zeros(T, 1)
    maxcut = maximum([model.rcut for model in model_matrix])
    cell_list = list_type(box, maxcut, N)
    system = Molecules(position, species, molecule, molecule_species, density, temperature, energy, model_matrix, d, N, box, e_locals, cell_list, bonds)
    build_cell_list!(system, system.cell_list)
    system.local_energy .= [compute_local_energy(system, i, cell_list) for i in particle_ids]
    system.energy[1] = sum(system.local_energy) / 2
    return system
end

function check_compute_energy_ij(system::Molecules, i, j, position_i, model_ij::Model, bonded)
    # Early return using && for short-circuit evaluation
    i == j && return zero(typeof(system.density))
    # If i != j, compute energy directly
    position_j = get_position(system, j)
    return compute_energy_ij(system, position_i, position_j, model_ij, bonded)
end

function compute_energy_ij(system::Molecules, position_i, position_j, model_ij::Model, bonded)
    energy_ij = zero(typeof(system.density))
    r2 = nearest_image_distance_squared(position_i, position_j, get_box(system))
    if r2 ≤ cutoff2(model_ij)
        energy_ij += potential(r2, model_ij)
    end
    if bonded
        energy_ij += bond_potential(r2, model_ij)
    end
    return energy_ij
end

function compute_energy_bonded_i(system::Molecules, i, position_i, bonds_i)
    energy_bonded_i = zero(typeof(system.density))
    for j in bonds_i
        model_ij = get_model(system, i, j)
        energy_bonded_i += compute_energy_ij(system, position_i, get_position(system, j), model_ij, true)
    end
    return energy_bonded_i
end

function compute_local_energy(system::Molecules, i, ::EmptyList)
    energy = zero(typeof(system.density))
    position_i = system.position[i]
    bonds_i = system.bonds[i]
    @inbounds for j in system.particle_ids
        bonded = j ∈ bonds_i
        model_ij = get_model(system, i, j)
        energy += check_compute_energy_ij(system, i, j, position_i, model_ij, bonded)
    end
    return energy
end

function destroy_particle!(system::Molecules, i, ::EmptyList)
    energy = zero(typeof(system.density))
    position_i = system.position[i]
    bonds_i = system.bonds[i]
    # Loop over particles
    @inbounds for j in system.particle_ids
        bonded = j ∈ bonds_i
        model_ij = get_model(system, i, j)
        energy_ij = check_compute_energy_ij(system, i, j, position_i, model_ij, bonded)
        energy += energy_ij
    end
    return energy
end

function create_particle!(system::Molecules, i, ::EmptyList)
    energy = zero(typeof(system.density))
    position_i = get_position(system, i)
    bonds_i = system.bonds[i]
    # Loop over particles
    @inbounds for j in system.particle_ids
        bonded = j ∈ bonds_i
        model_ij = get_model(system, i, j)
        energy_ij = check_compute_energy_ij(system, i, j, position_i, model_ij, bonded)
        energy += energy_ij
    end
    return energy
end

function compute_local_energy(system::Molecules, i, cell_list::LinkedList)
    energy_i = zero(typeof(system.density))
    position_i = system.position[i]
    mc = get_cell(position_i, cell_list)
    # Scan the neighbourhood of cell mc (including itself)
    bonds_i = system.bonds[i]
    energy_i += compute_energy_bonded_i(system, i, position_i, bonds_i)
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(cell_list, mc2)
        # Scan atoms in cell c2
        j = cell_list.head[c2]
        while (j != -1)
            if j ∉ bonds_i
                model_ij = get_model(system, i, j)
                energy_i += check_compute_energy_ij(system, i, j, position_i, model_ij, false)
            end
            j = cell_list.list[j]
        end
    end
    return energy_i
end


# With linked list
function destroy_particle!(system::Molecules, i, cell_list::LinkedList)
    energy_i = zero(typeof(system.density))
    # Get cell of particle i
    position_i = system.position[i]
    mc = get_cell(position_i, cell_list)
    # Scan the neighbourhood of cell mc (including itself)
    bonds_i = system.bonds[i]
    energy_i += compute_energy_bonded_i(system, i, position_i, bonds_i)
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(cell_list, mc2)
        # Scan atoms in cell c2
        j = cell_list.head[c2]
        while (j != -1)
            if j ∉ bonds_i
                model_ij = get_model(system, i, j)
                energy_ij = check_compute_energy_ij(system, i, j, position_i, model_ij, false)
                energy_i += energy_ij
            end
            j = cell_list.list[j]
        end
    end
    return energy_i
end

function create_particle!(system::Molecules, i, cell_list::LinkedList)
    energy_i = zero(typeof(system.density))
    # Get cell of particle i
    position_i = system.position[i]
    mc = get_cell(position_i, cell_list)
    # Scan the neighbourhood of cell mc (including itself)
    bonds_i = system.bonds[i]
    energy_i += compute_energy_bonded_i(system, i, position_i, bonds_i)
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(cell_list, mc2)
        # Scan atoms in cell c2
        j = cell_list.head[c2]
        while (j != -1)
            if j ∉ bonds_i
                model_ij = get_model(system, i, j)
                energy_ij = check_compute_energy_ij(system, i, j, position_i, model_ij, false)
                energy_i += energy_ij
            end
            j = cell_list.list[j]
        end
    end
    return energy_i
end


function compute_local_energy(system::Molecules, i, cell_list::CellList)
    energy_i = zero(typeof(system.density))
    position_i = get_position(system, i)
    mc = get_cell(position_i, cell_list)
    # Scan the neighbourhood of cell mc (including itself)
    bonds_i = system.bonds[i]
    energy_i += compute_energy_bonded_i(system, i, position_i, bonds_i)
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(cell_list, mc2)
        # Scan atoms in cell c2
        for j in cell_list.cells[c2]
            if j ∉ bonds_i
                model_ij = get_model(system, i, j)
                energy_i += check_compute_energy_ij(system, i, j, position_i, model_ij, false)
            end
        end
    end
    return energy_i
end

# With linked list
function destroy_particle!(system::Molecules, i, cell_list::CellList)
    # Get cell of particle i
    energy_i = zero(typeof(system.density))
    position_i = get_position(system, i)
    mc = get_cell(position_i, cell_list)
    bonds_i = system.bonds[i]
    energy_i += compute_energy_bonded_i(system, i, position_i, bonds_i)
    # Scan the neighbourhood of cell mc (including itself)
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(cell_list, mc2)
        # Scan atoms in cell c2
        neighbours = cell_list.cells[c2]
        @inbounds for j in neighbours
            if j ∉ bonds_i
                model_ij = get_model(system, i, j)
                energy_i += check_compute_energy_ij(system, i, j, position_i, model_ij, false)
            end
        end
    end
    return energy_i
end

function create_particle!(system::Molecules, i, cell_list::CellList)
    energy_i = zero(typeof(system.density))
    # Get cell of particle i
    position_i = get_position(system, i)
    mc = get_cell(position_i, cell_list)
    bonds_i = system.bonds[i]
    energy_i += compute_energy_bonded_i(system, i, position_i, bonds_i)
    # Scan the neighbourhood of cell mc (including itself)
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(cell_list, mc2)
        # Scan atoms in cell c2
        neighbours = cell_list.cells[c2]
        @inbounds for j in neighbours
            if j ∉ bonds_i
                model_ij = get_model(system, i, j)
                energy_i += check_compute_energy_ij(system, i, j, position_i, model_ij, false)
            end
        end
    end
    return energy_i
end
