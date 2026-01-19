"""
Atoms - type representing a system of non-bonded particles.

Description
`Atoms{D,V,C,H,T,SM} <: Particles` stores positions, species, thermodynamic
properties, and neighbour-list data needed for energy and neighbour queries.

Fields
- `position::Vector{SVector{D,T}}`: positions of particles.
- `species::V`: species/type identifier per particle (or per-site values).
- `density::T`, `temperature::T`, `energy::Vector{T}`: thermodynamic properties.
- `model_matrix::SM`: interaction models used for pairwise energies.
- `N::Int`, `d::Int`: number of particles and dimensionality.
- `box::SVector{D,T}`: periodic box vector.
- `neighbour_list::C`: neighbour list structure (CellList/LinkedList/EmptyList).
- `species_list::H`: optional auxiliary species list structure.
"""
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


"""
Create an `Atoms` system, initialize its neighbour list, and compute initial energy.

`System(position, species, density, temperature, model_matrix; list_type=EmptyList)`
constructs an `Atoms` instance, builds the neighbour list of type `list_type`,
and computes the initial total energy (stored in `system.energy[1]`).
"""
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

"""
Compute the pair energy between particles `i` and `j` for an `Atoms` system.

Returns zero if `i == j` or if the squared distance exceeds the model's cutoff.
Otherwise returns `potential(r2, model_ij)` where `r2` is the squared nearest-image distance.
"""
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


"""
Compute the energy of particle `i` by brute force (no neighbour list).

`compute_energy_particle(system, i, ::EmptyList)` sums interactions of particle `i`
with all other particles using `compute_energy_ij`.
"""
function compute_energy_particle(system::Atoms, i, neighbour_list::EmptyList)
    energy_i = zero(typeof(system.density))
    position_i = get_position(system, i)
    for j in get_neighbour_indices(system, neighbour_list, i)
        energy_i += compute_energy_ij(system, i, j, position_i)
    end
    return energy_i
end


"""
Compute the energy of particle `i` using a `CellList` neighbour list.

This restricts pair evaluations to particles in neighbouring cells of `i`.
"""
function compute_energy_particle(system::Atoms, i, neighbour_list::CellList)
    energy_i = zero(typeof(system.density))
    position_i = get_position(system, i)
    for j in get_neighbour_indices(system, neighbour_list, i)
        energy_i += compute_energy_ij(system, i, j, position_i)
    end
    return energy_i
end

"""
Compute the energy of particle `i` using a `LinkedList` neighbour list.

This variant iterates linked list heads for neighbouring cells and accumulates
pair energies computed with `compute_energy_ij`.
"""
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
