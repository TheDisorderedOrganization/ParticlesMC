struct Molecules{D,  VS<:AbstractVector, C<:NeighbourList, T<:AbstractFloat, SM<:AbstractArray} <: Particles
    position::Vector{SVector{D,T}}
    species::VS
    molecule::VS
    molecule_species::VS
    start_mol::Vector{Int}
    length_mol::Vector{Int}
    density::T
    temperature::T
    energy::Vector{T}
    model_matrix::SM
    d::Int
    N::Int
    Nmol::Int
    box::SVector{D,T}
    local_energy::Vector{T}
    cell_list::C
    bonds::Vector{Vector{Int}}
end


function System(position, species, molecule, density::T, temperature::T, model_matrix, bonds; molecule_species=nothing, list_type=EmptyList) where {T<:AbstractFloat}
    @assert length(position) == length(species)
    N = length(position)
    Nmol = length(unique(molecule))
    start_mol, length_mol = get_first_and_counts(molecule)
    molecule_species = something(molecule_species, ones(Int, N))
    d = length(Array(position)[1])
    box = @SVector fill(T((N / density)^(1 / d)), d)
    local_energy = zeros(T, N)
    energy = zeros(T, 1)
    maxcut = maximum([model.rcut for model in model_matrix])
    cell_list = list_type(box, maxcut, N)
    system = Molecules(position, species, molecule, molecule_species,  start_mol, length_mol, density, temperature, energy, model_matrix, d, N, Nmol,box, local_energy, cell_list, bonds)
    build_cell_list!(system, system.cell_list)
    system.local_energy .= [compute_local_energy(system, i, cell_list) for i in eachindex(position)]
    system.energy[1] = sum(system.local_energy) / 2
    return system
end

get_start_end_mol(system::Particles, i::Int) = system.start_mol[i], system.start_mol[i] + system.length_mol[i] - 1

function get_first_and_counts(vec::Vector{Int})
    firsts = Int[]
    counts = Int[]
    
    # Handle empty vector case
    isempty(vec) && return firsts, counts
    
    # Initialize with first element
    current = vec[1]
    push!(firsts, 1)
    count = 1
    
    # Scan through vector
    @inbounds for i in 2:length(vec)
        if vec[i] != current
            push!(counts, count)
            push!(firsts, i)
            current = vec[i]
            count = 1
        else
            count += 1
        end
    end
    # Add last count
    push!(counts, count)
    
    return firsts, counts
end

function check_nonbonded_compute_energy_ij(system::Molecules, i, j, position_i, bonds_i)
    # Early return using && for short-circuit evaluation
    i == j && return zero(system.density)
    j ∈ bonds_i && return zero(system.density)
    # If i != j, compute energy directly
    position_j = get_position(system, j)
    model_ij = get_model(system, i, j)
    return compute_energy_ij(system, position_i, position_j, model_ij, false)
end

function check_compute_energy_ij(system::Molecules, i, j, position_i, bonds_i)
    # Early return using && for short-circuit evaluation
    i == j && return zero(system.density)
    bonded = j ∈ bonds_i
    # If i != j, compute energy directly
    position_j = get_position(system, j)
    model_ij = get_model(system, i, j)
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
    @inbounds for j in bonds_i
        energy_bonded_i += compute_energy_ij(system, position_i, get_position(system, j), get_model(system, i, j), true)
    end
    return energy_bonded_i
end

function compute_local_energy(system::Molecules, i, ::EmptyList)
    energy = zero(typeof(system.density))
    position_i = system.position[i]
    bonds_i = system.bonds[i]
    @inbounds for j in eachindex(system)
        energy += check_compute_energy_ij(system, i, j, position_i, bonds_i)
    end
    return energy
end

function destroy_particle!(system::Molecules, i, ::EmptyList)
    energy = zero(typeof(system.density))
    position_i = system.position[i]
    bonds_i = system.bonds[i]
    # Loop over particles
    @inbounds for j in eachindex(system)
        energy += check_compute_energy_ij(system, i, j, position_i, bonds_i)
    end
    return energy
end

function create_particle!(system::Molecules, i, ::EmptyList)
    energy_i = zero(typeof(system.density))
    position_i = get_position(system, i)
    bonds_i = system.bonds[i]
    # Loop over particles
    @inbounds for j in eachindex(system)
        energy_i += check_compute_energy_ij(system, i, j, position_i, bonds_i)
    end
    return energy_i
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
            energy_i += check_nonbonded_compute_energy_ij(system, i, j, position_i, bonds_i)
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
    c = cell_index(cell_list, mc)
    # Scan the neighbourhood of cell mc (including itself)
    bonds_i = system.bonds[i]
    energy_i += compute_energy_bonded_i(system, i, position_i, bonds_i)
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(cell_list, mc2)
        j = cell_list.head[c2]
        while (j != -1)
            energy_i += check_nonbonded_compute_energy_ij(system, i, j, position_i, bonds_i)
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
            energy_i += check_nonbonded_compute_energy_ij(system, i, j, position_i, bonds_i)
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
            energy_i += check_nonbonded_compute_energy_ij(system, i, j, position_i, bonds_i)
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
            energy_i += check_nonbonded_compute_energy_ij(system, i, j, position_i, bonds_i)
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
            energy_i += check_nonbonded_compute_energy_ij(system, i, j, position_i, bonds_i)

        end
    end
    return energy_i
end
