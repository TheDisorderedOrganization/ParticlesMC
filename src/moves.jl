using ComponentArrays

"""
Execute an `Action` via the Arianna API and update the system energy.

`Arianna.perform_action!(system, action)` calls `perform_action!` to obtain the
pre- and post-action energies `(e₁, e₂)`, sets `action.δe` and updates
`system.energy[1]` accordingly (ignoring updates when energies are Inf).
It returns the tuple `(e₁, e₂)`.
"""
function Arianna.perform_action!(system::Particles, action::Action)
    e₁, e₂ = perform_action!(system, action)
    if isinf(e₁) || isinf(e₂)
        action.δe = zero(typeof(system.energy[1]))
    else
        action.δe = e₂ - e₁
        system.energy[1] += action.δe
    end
    return e₁, e₂
end

###############################################################################
# SIMPLE DISPLACEMENT

"""
    mutable struct Displacement{T<:AbstractArray} <: Action

A struct representing a displacement action, where particle i is moved by specified amounts `δ`.

# Fields
- `i::Int`: Indices of the particles or elements in `system` to be displaced.
- `δ::T`: Displacement values for each corresponding index in `is`.
"""
mutable struct Displacement{T<:AbstractArray, F<:AbstractFloat} <: Action
    i::Int
    δ::T
    δe::F
end

"""
Apply the displacement contained in `action` to the particle(s) in `system`.

`update_position!(system, action::Displacement)` updates `system.position[action.i]`
by adding `action.δ` in-place.
"""
function update_position!(system::Particles, action::Displacement)
    @inbounds system.position[action.i] = system.position[action.i] + action.δ
end

"""
Perform a displacement action: compute pre-move energy, apply the displacement,
update the neighbour list if the particle changed cell, then compute post-move
energy and store the energy change in `action.δe`.

Returns the tuple `(e₁, e₂)` of energies before and after the move.
"""
function perform_action!(system::Particles, action::Displacement)
    neighbour_list = get_neighbour_list(system)
    e₁ = compute_energy_particle(system, action.i, neighbour_list)
    update_position!(system, action)
    c, c2 = old_new_cell(system, action.i, neighbour_list)
    if c != c2
        update_neighbour_list!(action.i, c, c2, neighbour_list)
    end
    e₂ = compute_energy_particle(system, action.i, neighbour_list)
    action.δe = e₂ - e₁
    return e₁, e₂
end

"""
Revert a previously applied `Displacement` action.

This function applies `action.δ` (which should typically be the inverse
displacement), restores the system energy by subtracting `action.δe`, and
updates the neighbour list if the particle changed cell.
"""
function Arianna.revert_action!(system::Particles, action::Displacement)
    update_position!(system, action)
    neighbour_list = get_neighbour_list(system)
    system.energy[1] -= action.δe
    c, c2 = old_new_cell(system, action.i, neighbour_list)
    if c != c2
        update_neighbour_list!(action.i, c, c2, neighbour_list)
    end
end

"""
Invert a `Displacement` action in-place by negating its displacement vector.

This is used to create the revert operation for a displacement move.
"""
function Arianna.invert_action!(action::Displacement, ::Particles)
    action.δ = -action.δ
end

"""
Reward function for `Displacement` under a policy-guided scheme.

Returns the squared norm of the displacement (can be used as a positive reward
for exploration).
"""
function Arianna.PolicyGuided.reward(action::Displacement, ::Particles)
    return norm(action.δ)^ 2
end

"""
Simple isotropic Gaussian proposal policy for displacement moves.
"""
struct SimpleGaussian <: Policy end

"""
Compute the log proposal density of a `Displacement` under `SimpleGaussian`.
"""
function Arianna.log_proposal_density(action::Displacement, ::SimpleGaussian, parameters, system::Particles)
    return -norm(action.δ)^2 / (2parameters.σ^2) - system.d * log(2π * parameters.σ^2) / 2
end

"""
Sample a `Displacement` action under `SimpleGaussian`.

Chooses a random particle index and draws an isotropic Gaussian displacement of
dimension `system.d` with standard deviation `parameters.σ`.
"""
function Arianna.sample_action!(action::Displacement, ::SimpleGaussian, parameters, system::Particles, rng)
    action.i = rand(rng, 1:length(system))
    action.δ = rand(rng, Normal(0, parameters.σ), system.d)
end

###############################################################################
# DISCRETESWAP

"""
DiscreteSwap - action representing swapping species between two particles.

Fields
- `i::Int`, `j::Int`: particle indices involved in the swap.
- `species::Tuple{Int,Int}`: pair of species identifiers being swapped.
- `particles_per_species::Tuple{Int,Int}`: population counts used for proposal densities.
- `δe::Float64`: energy difference produced by the swap.
"""
mutable struct DiscreteSwap <: Action
    i::Int
    j::Int
    species::Tuple{Int,Int}
    particles_per_species::Tuple{Int,Int}
    δe::Float64
end

function DiscreteSwap(species::Vector{Int}, system::Atoms)
    particles_per_species_1 = count(i -> (i == species[1]), system.species)
    particles_per_species_2 = count(i -> (i == species[2]), system.species)
    return DiscreteSwap(1, 1, (species[1], species[2]), (particles_per_species_1, particles_per_species_2), 0.0)
end

"""
Swap the species identifiers of particles `i` and `j` and return the total
pre- and post-swap energies for those two particles.

`swap_particle_species!(system, spi, i, spj, j)` computes energies before the
swap, performs the species exchange, then computes energies after and returns
`(E_before, E_after)` for the pair.
"""
function swap_particle_species!(system::Particles, spi, i, spj, j)
    e₁ᵢ = compute_energy_particle(system, i, system.neighbour_list)
    e₁ⱼ = compute_energy_particle(system, j, system.neighbour_list)
    system.species[i] = spj
    system.species[j] = spi
    e₂ᵢ = compute_energy_particle(system, i, system.neighbour_list)
    e₂ⱼ = compute_energy_particle(system, j, system.neighbour_list)
    return e₁ᵢ + e₁ⱼ, e₂ᵢ + e₂ⱼ
end

"""
Update auxiliary species list metadata after swapping particles `i` and `j`.

`update_species_list!(species_list, swap_species, i, j)` updates the internal
sp_ids and sp_heads structures to reflect the swap between indices `i` and `j`.
"""
function update_species_list!(species_list, swap_species, i, j)
    species_list.sp_ids[swap_species[1]][species_list.sp_heads[i]] = j
    species_list.sp_ids[swap_species[2]][species_list.sp_heads[j]] = i
    species_list.sp_heads[i], species_list.sp_heads[j] = species_list.sp_heads[j], species_list.sp_heads[i]
end

"""
Execute a `DiscreteSwap` action: swap species of particles `i` and `j`, store
`action.δe`, and update the species list auxiliary structure.

Returns `(E_before, E_after)` for the swapped particles.
"""
function perform_action!(system::Particles, action::DiscreteSwap)
    i, j = action.i, action.j
    spi, spj = get_species(system, i), get_species(system, j)
    e₁, e₂ = swap_particle_species!(system, spi, i, spj, j)
    action.δe = e₂ - e₁
    update_species_list!(system.species_list, action.species, i, j)
    return e₁, e₂
end

"""
Revert a `DiscreteSwap` action by swapping species back and restoring energy.

This function also updates the species list metadata to undo the swap.
"""
function Arianna.revert_action!(system::Particles, action::DiscreteSwap)
    i, j = action.i, action.j
    spi, spj = get_species(system, i), get_species(system, j)
    system.species[j], system.species[i] = spi, spj
    update_species_list!(system.species_list, action.species, i, j)
    system.energy[1] -= action.δe
end

"""
Invert a `DiscreteSwap` action by swapping the stored particle indices.
"""
function Arianna.invert_action!(action::DiscreteSwap, ::Particles)
    action.i, action.j = action.j, action.i
end

"""
Policy-guided reward function for `DiscreteSwap` (constant reward).
"""
function Arianna.PolicyGuided.reward(::DiscreteSwap, system::Particles)
    return one(typeof(system.temperature))
end

"""
DoubleUniform - proposal policy selecting one particle uniformly from each of two species lists.
"""
struct DoubleUniform <: Policy end

"""
Log proposal density for `DiscreteSwap` under `DoubleUniform` sampling.
"""
function Arianna.log_proposal_density(action::DiscreteSwap, ::DoubleUniform, parameters, system::Particles)
    return -log(action.particles_per_species[1] * action.particles_per_species[2])
end

"""
Sample indices for `DiscreteSwap` uniformly from the two species lists.
"""
function Arianna.sample_action!(action::DiscreteSwap, ::DoubleUniform, parameters, system::Particles, rng)
    action.i = rand(rng, system.species_list.sp_ids[action.species[1]])
    action.j = rand(rng, system.species_list.sp_ids[action.species[2]])
end

"""
EnergyBias - proposal policy biasing selection by particle energies.
"""
struct EnergyBias <: Policy end

"""
Log proposal density for `DiscreteSwap` under the `EnergyBias` policy.

This function computes the log-probability of selecting the pair `(i,j)` given
energy-weighted exponentials controlled by `parameters.θ₁` and `parameters.θ₂`.
"""
function Arianna.log_proposal_density(action::DiscreteSwap, ::EnergyBias, parameters, system::Particles)
    energy_particle_i = compute_energy_particle(system, action.i)
    energy_particle_j = compute_energy_particle(system, action.j)
    numerator = parameters.θ₁ * energy_particle_i+ parameters.θ₂ * energy_particle_j
    energy_particles_1 = compute_energy_particle(system, system.species_list.sp_ids[action.species[1]])
    energy_particles_2 = compute_energy_particle(system, system.species_list.sp_ids[action.species[2]])
    log_sum_exp_1 = log(sum(exp.(parameters.θ₁ .* energy_particles_1)))
    log_sum_exp_2 = log(sum(exp.(parameters.θ₂ .* energy_particles_2)))
    return numerator - log_sum_exp_1 - log_sum_exp_2
end

"""
Sample a `DiscreteSwap` pair under the `EnergyBias` policy by importance sampling
based on particle energies.
"""
function Arianna.sample_action!(action::DiscreteSwap, ::EnergyBias, parameters, system::Particles, rng)
    energy_particles_1 = compute_energy_particle(system, system.species_list.sp_ids[action.species[1]])
    energy_particles_2 = compute_energy_particle(system, system.species_list.sp_ids[action.species[2]])
    w1s = exp.(parameters.θ₁ .* energy_particles_1)
    w2s = exp.(parameters.θ₂ .* energy_particles_2)
    w1s .= w1s ./ sum(w1s)
    w2s .= w2s ./ sum(w2s)
    id1 = rand(rng, Categorical(w1s))
    id2 = rand(rng, Categorical(w2s))
    action.i = system.species_list.sp_ids[action.species[1]][id1]
    action.j = system.species_list.sp_ids[action.species[2]][id2]
end


###############################################################################
"""
MoleculeFlip - action representing swapping species between two sites within the same molecule.

Fields
- `i::Int`, `j::Int`: particle indices inside the same molecule to be swapped.
- `δe::F`: energy change due to the swap.
"""
mutable struct MoleculeFlip{F<:AbstractFloat} <: Action
    i::Int
    j::Int
    δe::F
end

"""
Perform a `MoleculeFlip` action by swapping species between two sites within a molecule.
Returns the pair energy before and after the swap `(e₁, e₂)` and stores `δe`.
"""
function perform_action!(system::Particles, action::MoleculeFlip)
    i, j = action.i, action.j
    spi, spj = system.species[i], system.species[j]
    e₁, e₂ = swap_particle_species!(system, spi, i, spj, j)
    action.δe = e₂ - e₁
    return e₁, e₂
end

"""
Invert a `MoleculeFlip` action by swapping stored indices.
"""
function Arianna.invert_action!(action::MoleculeFlip, ::Molecules)
    action.i, action.j = action.j, action.i
end

"""
Revert a `MoleculeFlip` action by swapping the species back and restoring energy.
"""
function Arianna.revert_action!(system::Particles, action::MoleculeFlip)
    i, j = action.i, action.j
    spi, spj = get_species(system, i), get_species(system, j)
    system.species[j], system.species[i] = spi, spj
    system.energy[1] -= action.δe
end

"""
Reward function for `MoleculeFlip` (currently constant).
"""
function reward(::MoleculeFlip, system::Particles)
    return one(typeof(system.temperature))
end

"""
Log proposal density for `MoleculeFlip` under `DoubleUniform` sampling (uniform between 2 choices).
"""
function Arianna.log_proposal_density(action::MoleculeFlip, ::DoubleUniform, parameters, system::Particles)
    return -log(2)
end

"""
Sample a `MoleculeFlip` action by selecting a molecule uniformly and picking two distinct
sites in that molecule with different species.
"""
function Arianna.sample_action!(action::MoleculeFlip, ::DoubleUniform, parameters, system::Particles, rng)
    molecule_i = rand(rng, DiscreteUniform(1, system.Nmol))
    start_mol, end_mol = get_start_end_mol(system, molecule_i)
    i, j = sample(rng, start_mol:end_mol, 2; replace=false)
    while system.species[i] == system.species[j]
        i, j = sample(rng, start_mol:end_mol, 2; replace=false)
    end
    action.i, action.j = i, j
end

###############################################################################
