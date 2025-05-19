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
    neighbour_list::C
    bonds::Vector{Vector{Int}}
end


struct Bonded end
struct NonBonded end

function System(position, species, molecule, density::T, temperature::T, model_matrix, bonds; molecule_species=nothing, list_type=EmptyList) where {T<:AbstractFloat}
    @assert length(position) == length(species)
    N = length(position)
    Nmol = length(unique(molecule))
    start_mol, length_mol = get_first_and_counts(molecule)
    molecule_species = something(molecule_species, ones(Int, N))
    d = length(Array(position)[1])
    box = @SVector fill(T((N / density)^(1 / d)), d)
    energy = zeros(T, 1)
    maxcut = maximum([model.rcut for model in model_matrix])
    neighbour_list = list_type(box, maxcut, N)
    system = Molecules(position, species, molecule, molecule_species,  start_mol, length_mol, density, temperature, energy, model_matrix, d, N, Nmol,box, neighbour_list, bonds)
    build_neighbour_list!(system)
    local_energy = [compute_energy_particle(system, i, neighbour_list) for i in eachindex(position)]
    energy = sum(local_energy) / 2
    if isinf(energy) || isnan(energy)
        error("Initial configuration has infinite or NaN energy.")
    end
    system.energy[1] = energy
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

function check_compute_energy_ij(system::Molecules, i, j, position_i, bonds_i)
    # Early return using && for short-circuit evaluation
    i == j && return zero(system.density)
    isbonded = j ∈ bonds_i ? Bonded() : NonBonded()
    # If i != j, compute energy directly
    position_j = get_position(system, j)
    model_ij = get_model(system, i, j)
    return compute_energy_ij(system, position_i, position_j, model_ij, isbonded)
end

function check_nonbonded_compute_energy_ij(system::Molecules, i, j, position_i, bonds_i)
    # Early return using && for short-circuit evaluation
    i == j && return zero(system.density)
    j ∈ bonds_i && return zero(system.density)
    # If i != j, compute energy directly
    position_j = get_position(system, j)
    model_ij = get_model(system, i, j)
    return compute_energy_ij(system, position_i, position_j, model_ij, NonBonded())
end

function compute_energy_bonded_i(system::Molecules, i, position_i, bonds_i)
    energy_bonded_i = zero(typeof(system.density))
    @inbounds for j in bonds_i
        energy_bonded_i += compute_energy_ij(system, position_i, get_position(system, j), get_model(system, i, j), Bonded())
    end
    return energy_bonded_i
end

function compute_energy_ij(system::Molecules, position_i, position_j, model_ij::Model, ::Bonded)
    box = get_box(system)
    r2 = nearest_image_distance_squared(position_i, position_j, box)
    energy_ij = bond_potential(r2, model_ij)
    cutoff2_val = cutoff2(model_ij)
    if r2 ≤ cutoff2_val
        energy_ij += potential(r2, model_ij)
    end

    return energy_ij
end
function compute_energy_ij(system::Molecules, position_i, position_j, model_ij::Model, ::NonBonded)
    r2 = nearest_image_distance_squared(position_i, position_j, get_box(system))
    cutoff2_val = cutoff2(model_ij)
    return r2 > cutoff2_val ? zero(typeof(system.density)) : potential(r2, model_ij)
end

function compute_energy_particle(system::Molecules, i, ::EmptyList)
    energy = zero(typeof(system.density))
    position_i = system.position[i]
    bonds_i = system.bonds[i]
    @inbounds for j in eachindex(system)
        energy += check_compute_energy_ij(system, i, j, position_i, bonds_i)
    end
    return energy
end

# With linked list
function compute_energy_particle(system::Molecules, i, neighbour_list::LinkedList)
    energy_i = zero(typeof(system.density))
    # Get cell of particle i
    position_i = system.position[i]
    c = get_cell_index(i, neighbour_list)
    cells = neighbour_list.neighbour_cells[c]
    # Scan the neighbourhood of cell mc (including itself)
    bonds_i = system.bonds[i]
    energy_i += compute_energy_bonded_i(system, i, position_i, bonds_i)
    @inbounds for c2 in cells
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        j = neighbour_list.head[c2]
        while (j != -1)
            energy_i += check_nonbonded_compute_energy_ij(system, i, j, position_i, bonds_i)
            j = neighbour_list.list[j]
        end
    end
    return energy_i
end

function compute_energy_particle(system::Molecules, i, neighbour_list::CellList)
    energy_i = zero(typeof(system.density))
    position_i = get_position(system, i)
    c = get_cell_index(i, neighbour_list)
    neighbour_cells = neighbour_list.neighbour_cells[c]
    # Scan the neighbourhood of cell mc (including itself)
    bonds_i = system.bonds[i]
    energy_i += compute_energy_bonded_i(system, i, position_i, bonds_i)
    @inbounds for c2 in neighbour_cells
        # Scan atoms in cell c2
        neighbours = neighbour_list.cells[c2]
        @inbounds for j in neighbours
            energy_i += check_nonbonded_compute_energy_ij(system, i, j, position_i, bonds_i)
        end
    end
    return energy_i
end

function compute_chain_correlation(system::Molecules)
    chain_length = system.length_mol[1]
    @assert all(x -> x == chain_length, system.length_mol) "All chains must have the same length"
    @assert chain_length > 1 "Chains must have at least two particles"
    polymer_array = zeros(system.Nmol, chain_length)
    for i in 1:system.Nmol
        start = system.start_mol[i]
        polymer_array[i, :] = system.species[start:start+chain_length-1]
    end
    polymer_array[polymer_array .== 2] .= -1
    correlation_array = Float64[]
    for i in 1:chain_length-1
        for j in i+1:chain_length
            cross = sum(polymer_array[:, i] .* polymer_array[:, j]) / system.Nmol
            push!(correlation_array, cross)
        end
    end
    return sum(correlation_array.^2)
end