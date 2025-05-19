using ComponentArrays

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

function update_position!(system::Particles, action::Displacement)
    @inbounds system.position[action.i] = system.position[action.i] + action.δ
end

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

function Arianna.revert_action!(system::Particles, action::Displacement)
    update_position!(system, action)
    neighbour_list = get_neighbour_list(system)
    system.energy[1] -= action.δe
    c, c2 = old_new_cell(system, action.i, neighbour_list)
    if c != c2
        update_neighbour_list!(action.i, c, c2, neighbour_list)
    end
end

function Arianna.invert_action!(action::Displacement, ::Particles)
    action.δ = -action.δ
end

function Arianna.PolicyGuided.reward(action::Displacement, ::Particles)
    return norm(action.δ)^ 2
end

struct SimpleGaussian <: Policy end

function Arianna.log_proposal_density(action::Displacement, ::SimpleGaussian, parameters, system::Particles)
    return -norm(action.δ)^2 / (2parameters.σ^2) - system.d * log(2π * parameters.σ^2) / 2
end

function Arianna.sample_action!(action::Displacement, ::SimpleGaussian, parameters, system::Particles, rng)
    action.i = rand(rng, 1:length(system))
    action.δ = rand(rng, Normal(0, parameters.σ), system.d)
end

###############################################################################
# DISCRETESWAP

mutable struct DiscreteSwap <: Action
    i::Int
    j::Int
    species::Tuple{Int,Int}
    particles_per_species::Tuple{Int,Int}
    δe::Float64
end

function swap_particle_species!(system::Particles, spi, i, spj, j)
    e₁ᵢ = compute_energy_particle(system, i, system.neighbour_list)
    e₁ⱼ = compute_energy_particle(system, j, system.neighbour_list)
    system.species[i] = spj
    system.species[j] = spi
    e₂ᵢ = compute_energy_particle(system, i, system.neighbour_list)
    e₂ⱼ = compute_energy_particle(system, j, system.neighbour_list)
    return e₁ᵢ + e₁ⱼ, e₂ᵢ + e₂ⱼ
end

function update_species_list!(species_list, swap_species, i, j)
    species_list.sp_ids[swap_species[1]][species_list.sp_heads[i]] = j
    species_list.sp_ids[swap_species[2]][species_list.sp_heads[j]] = i
    species_list.sp_heads[i], species_list.sp_heads[j] = species_list.sp_heads[j], species_list.sp_heads[i]
end

function perform_action!(system::Particles, action::DiscreteSwap)
    i, j = action.i, action.j
    spi, spj = get_species(system, i), get_species(system, j)
    e₁, e₂ = swap_particle_species!(system, spi, i, spj, j)
    action.δe = e₂ - e₁
    update_species_list!(system.species_list, action.species, i, j)
    return e₁, e₂
end

function Arianna.revert_action!(system::Particles, action::DiscreteSwap)
    i, j = action.i, action.j
    spi, spj = get_species(system, i), get_species(system, j)
    system.species[j], system.species[i] = spi, spj
    update_species_list!(system.species_list, action.species, i, j)
    system.energy[1] -= action.δe
end

function Arianna.invert_action!(action::DiscreteSwap, ::Particles)
    action.i, action.j = action.j, action.i
end

function Arianna.PolicyGuided.reward(::DiscreteSwap, system::Particles)
    return one(typeof(system.temperature))
end

struct DoubleUniform <: Policy end

function Arianna.log_proposal_density(action::DiscreteSwap, ::DoubleUniform, parameters, system::Particles)
    return -log(action.particles_per_species[1] * action.particles_per_species[2])
end

function Arianna.sample_action!(action::DiscreteSwap, ::DoubleUniform, parameters, system::Particles, rng)
    action.i = rand(rng, system.species_list.sp_ids[action.species[1]])
    action.j = rand(rng, system.species_list.sp_ids[action.species[2]])
end

struct EnergyBias <: Policy end

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
mutable struct MoleculeFlip{F<:AbstractFloat} <: Action
    i::Int
    j::Int
    δe::F
end

function perform_action!(system::Particles, action::MoleculeFlip)
    i, j = action.i, action.j
    spi, spj = system.species[i], system.species[j]
    e₁, e₂ = swap_particle_species!(system, spi, i, spj, j)
    action.δe = e₂ - e₁
    return e₁, e₂
end

function Arianna.invert_action!(action::MoleculeFlip, ::Molecules)
    action.i, action.j = action.j, action.i
end

function Arianna.revert_action!(system::Particles, action::MoleculeFlip)
    i, j = action.i, action.j
    spi, spj = get_species(system, i), get_species(system, j)
    system.species[j], system.species[i] = spi, spj
    system.energy[1] -= action.δe
end

function reward(::MoleculeFlip, system::Particles)
    return one(typeof(system.temperature))
end

function Arianna.log_proposal_density(action::MoleculeFlip, ::DoubleUniform, parameters, system::Particles)
    return -log(2)
end

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