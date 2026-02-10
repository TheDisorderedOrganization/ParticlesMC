"""
Molecules - type managing a collection of molecules (chains) of particles.

Description
`Molecules{D,VS,C,T,SM} <: Particles` represents a particle system organized into
molecules (or chains). It stores particle positions, molecule membership,
molecule lengths and start indices, thermodynamic properties, interaction
models, and the neighbour list structure used for energy and neighbour queries.

Fields
- `position::Vector{SVector{D,T}}`: positions of particles.
- `species::VS`: species/type of each particle.
- `molecule::VS`: molecule identifier for each particle.
- `molecule_species::VS`: species identifier for each molecule.
- `start_mol::Vector{Int}`: starting site index of each molecule.
- `length_mol::Vector{Int}`: length (number of sites) of each molecule.
- `density::T`, `temperature::T`, `energy::Vector{T}`: thermodynamic properties.
- `model_matrix::SM`: interaction model matrix.
- `d::Int`, `N::Int`, `Nmol::Int`: dimension, number of particles, number of molecules.
- `box::SVector{D,T}`: periodic box vector.
- `neighbour_list::C`: neighbour list structure (CellList/LinkedList/EmptyList).
- `bonds::Vector{Vector{Int}}`: lists of bonded neighbours for each particle.
"""
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

"""
Marker indicating that a pair of particles is bonded.

Used to select the bonded contribution (bond potentials) in energy computations.
"""
struct Bonded end

"""
Marker indicating that a pair of particles is non-bonded.

Used to select the non-bonded pair potential contribution in energy computations.
"""
struct NonBonded end

"""
Create a `Molecules` system and initialize neighbour list and energy.

`System(position, species, molecule, density, temperature, model_matrix, bonds; molecule_species=nothing, list_type=EmptyList)` builds
and returns a `Molecules` instance. It computes `start_mol`/`length_mol`,
constructs the neighbour list of type `list_type`, builds it, and computes the
initial total energy (with a check for Inf/NaN).

Arguments
- `position`: vector of particle positions.
- `species`: species identifier per particle.
- `molecule`: molecule id per particle.
- `density`, `temperature`: thermodynamic scalars.
- `model_matrix`: interaction models used to compute energies.
- `bonds`: bonded neighbour lists for each particle.

Returns
- `Molecules` instance with neighbour list built and `energy[1]` set.
"""
function System(position, species, molecule, density::T, temperature::T, model_matrix, bonds; molecule_species=nothing, list_type=EmptyList, list_parameters=nothing) where {T<:AbstractFloat}
    @assert length(position) == length(species)
    N = length(position)
    Nmol = length(unique(molecule))
    start_mol, length_mol = get_first_and_counts(molecule)
    molecule_species = something(molecule_species, ones(Int, N))
    d = length(Array(position)[1])
    box = @SVector fill(T((N / density)^(1 / d)), d)
    energy = zeros(T, 1)
    maxcut = maximum([model.rcut for model in model_matrix])
    neighbour_list = list_type(box, maxcut, N; list_parameters=list_parameters)
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

"""
Return the start and end indices of molecule `i` in `system`.

`get_start_end_mol(system, i)` returns a tuple `(start, end)` for the i-th molecule,
where `start` is the first particle index and `end` is the last particle index.
"""
get_start_end_mol(system::Particles, i::Int) = system.start_mol[i], system.start_mol[i] + system.length_mol[i] - 1

"""
Return arrays of first indices and counts for consecutive blocks in `vec`.

`get_first_and_counts(vec)` scans `vec` and returns a tuple `(firsts, counts)`, where
`firsts` are the starting indices of each run and `counts` the corresponding lengths.
"""
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

"""
Check and compute the pair energy between particles `i` and `j`.

`check_compute_energy_ij(system, i, j, position_i, bonds_i)` returns the pair energy
between `i` and `j`. If `i == j` it returns zero. It determines whether the pair
is bonded and dispatches to `compute_energy_ij` with the appropriate marker
(`Bonded` or `NonBonded`).
"""
function check_compute_energy_ij(system::Molecules, i, j, position_i, bonds_i)
    # Early return using && for short-circuit evaluation
    i == j && return zero(system.density)
    isbonded = j ∈ bonds_i ? Bonded() : NonBonded()
    # If i != j, compute energy directly
    position_j = get_position(system, j)
    model_ij = get_model(system, i, j)
    return compute_energy_ij(system, position_i, position_j, model_ij, isbonded)
end
"""
`check_nonbonded_compute_energy_ij` returns zero if `i == j` or if `j` is bonded to `i`.
Otherwise it computes the non-bonded pair energy by dispatching to
`compute_energy_ij(..., NonBonded())`.
"""
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
    return bond_potential(r2, model_ij)
end

"""
Compute non-bonded pair energy using pair potential with cutoff.

This method is chosen when the pair marker is `NonBonded`. It returns zero if
the squared distance exceeds the model's cutoff squared, otherwise returns the
pair potential `potential(r2, model_ij)`.
"""
function compute_energy_ij(system::Molecules, position_i, position_j, model_ij::Model, ::NonBonded)
    r2 = nearest_image_distance_squared(position_i, position_j, get_box(system))
    cutoff2_val = cutoff2(model_ij)
    return r2 > cutoff2_val ? zero(typeof(system.density)) : potential(r2, model_ij)
end

"""
Compute particle energy using a the provided neighbour list.

Non-bonded pair energy evaluations are restricted to particles in the neighbour list;
bonded contributions are added explicitly.
"""
function compute_energy_particle(system::Molecules, i, neighbour_list::NeighbourList)
    position_i = system.position[i]
    bonds_i = system.bonds[i]

    energy_i = compute_energy_bonded_i(system, i, position_i, bonds_i)
    for j in neighbour_list(system, i)
        energy_i += check_nonbonded_compute_energy_ij(system, i, j, position_i, bonds_i)
    end
    return energy_i
end

"""
Compute a measure of chain correlation for monodisperse chains.

`compute_chain_correlation(system)` assumes all chains have the same length and
computes the squared sum of pairwise cross-correlations between different chain
sites over all chains (useful as an order parameter).
"""
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

Arianna.@callback function chain_correlation(system::Molecules)
    return compute_chain_correlation(system)
end
