using ComponentArrays


function cache_update!(system::Particles, action::Action)
    index_set = BitSet()
    @inbounds for elt in system.cache
        i = elt[1]
        if !(i in index_set)
            push!(index_set, i)
            system.local_energy[i] = elt[2] 
        end
    end
end

function build_cache!(system::Particles, i, ::EmptyList)
    return nothing
end

function build_cache!(system::Particles, i, cell_list::LinkedList)
    mc = get_cell(system.position[i], cell_list.cell)
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(mc2, cell_list.ncells)
        # Scan atoms in cell c2
        j = cell_list.head[c2]
        while (j != -1)
            push!(system.cache, (j, system.local_energy[j]))
            j = cell_list.list[j]
        end
    end
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


mutable struct Displacement{T<:AbstractArray} <: Action
    i::Int
    δ::T
end


function MonteCarlo.perform_action!(system::Particles, action::Displacement)
    empty!(system.cache)
    e₁ = destroy_particle!(system, action.i, system.cell_list)
    system.position[action.i] = system.position[action.i] + action.δ
    c, c2 = old_new_cell(system, action.i, system.cell_list)
    if c != c2
        update_cell_list!(system, action.i, c, c2, system.cell_list)
        build_cache!(system, action.i, system.cell_list)
    end
    e₂ = create_particle!(system, action.i, system.cell_list)
    return e₁, e₂
end

function MonteCarlo.perform_action_cached!(system::Particles, action::Displacement)
    system.position[action.i] = system.position[action.i] + action.δ
    c, c2 = old_new_cell(system, action.i, system.cell_list)
    if c != c2
        update_cell_list!(system, action.i, c, c2, system.cell_list)
    end
    cache_update!(system, action)
end

function MonteCarlo.invert_action!(action::Displacement, ::Particles)
    action.δ = -action.δ
    return nothing
end

function MonteCarlo.PolicyGuided.reward(action::Displacement, ::Particles)
    return norm(action.δ)^ 2
end

struct SimpleGaussian <: Policy end

function MonteCarlo.log_proposal_density(action::Displacement, ::SimpleGaussian, parameters, system::Particles)
    return -norm(action.δ)^2 / (2parameters.σ^2) - system.d * log(2π * parameters.σ^2) / 2
end

function MonteCarlo.sample_action!(action::Displacement, ::SimpleGaussian, parameters, system, rng)
    action.i = rand(rng, DiscreteUniform(1, system.N))
    action.δ = map(x -> rand(rng, Normal(x, parameters.σ)), zero(system.box))
    return nothing
end

###############################################################################
# DISCRETESWAP

mutable struct DiscreteSwap <: Action
    i::Int
    j::Int
    species::Tuple{Int,Int}
    particles_per_species::Tuple{Int,Int}
end

function swap_particle_species!(system::Particles, spi, i, spj, j)
    e₁ᵢ = destroy_particle!(system, i, system.cell_list)
    system.species[i] = spj
    e₂ᵢ = create_particle!(system, i, system.cell_list)
    e₁ⱼ = destroy_particle!(system, j, system.cell_list)
    system.species[j] = spi
    e₂ⱼ = create_particle!(system, j, system.cell_list)
    return e₁ᵢ + e₁ⱼ, e₂ᵢ + e₂ⱼ
end

function update_species_list!(species_list, swap_species, i, j)
    species_list.sp_ids[swap_species[1]][species_list.sp_heads[i]] = j
    species_list.sp_ids[swap_species[2]][species_list.sp_heads[j]] = i
    species_list.sp_heads[i], species_list.sp_heads[j] = species_list.sp_heads[j], species_list.sp_heads[i]
end

function MonteCarlo.perform_action!(system::Particles, action::DiscreteSwap)
    empty!(system.cache)
    i, j = action.i, action.j
    spi, spj = system.species[i], system.species[j]
    e₁, e₂ = swap_particle_species!(system, spi, i, spj, j)
    update_species_list!(system.species_list, action.species, i, j)
    return e₁, e₂
end

function MonteCarlo.perform_action_cached!(system::Particles, action::DiscreteSwap)
    i, j = action.i, action.j
    spi, spj = system.species[i], system.species[j]
    system.species[j], system.species[i] = spi, spj
    update_species_list!(system.species_list, action.species, i, j)
    cache_update!(system, action)
end

function MonteCarlo.invert_action!(action::DiscreteSwap, ::Particles)
    action.i, action.j = action.j, action.i
    return nothing
end

function MonteCarlo.PolicyGuided.reward(::DiscreteSwap, system::Particles)
    return one(typeof(system.temperature))
end

struct DoubleUniform <: Policy end

function MonteCarlo.log_proposal_density(action::DiscreteSwap, ::DoubleUniform, parameters, system::Particles)
    return -log(action.particles_per_species[1] * action.particles_per_species[2])
end

function MonteCarlo.sample_action!(action::DiscreteSwap, ::DoubleUniform, parameters, system::Particles, rng)
    action.i = rand(rng, system.species_list.sp_ids[action.species[1]])
    action.j = rand(rng, system.species_list.sp_ids[action.species[2]])
    return nothing
end

struct EnergyBias <: Policy end

function MonteCarlo.log_proposal_density(action::DiscreteSwap, ::EnergyBias, parameters, system::Particles)
    numerator = parameters.θ₁ * system.local_energy[action.i] + parameters.θ₂ * system.local_energy[action.j]
    log_sum_exp_1 = log(sum(exp.(parameters.θ₁ .* system.local_energy[system.species_list.sp_ids[action.species[1]]])))
    log_sum_exp_2 = log(sum(exp.(parameters.θ₂ .* system.local_energy[system.species_list.sp_ids[action.species[2]]])))
    return numerator - log_sum_exp_1 - log_sum_exp_2
end

function MonteCarlo.sample_action!(action::DiscreteSwap, ::EnergyBias, parameters, system::Particles, rng)
    w1s = exp.(parameters.θ₁ .* system.local_energy[system.species_list.sp_ids[action.species[1]]])
    w2s = exp.(parameters.θ₂ .* system.local_energy[system.species_list.sp_ids[action.species[2]]])
    w1s .= w1s ./ sum(w1s)
    w2s .= w2s ./ sum(w2s)
    id1 = rand(rng, Categorical(w1s))
    id2 = rand(rng, Categorical(w2s))
    action.i = system.species_list.sp_ids[action.species[1]][id1]
    action.j = system.species_list.sp_ids[action.species[2]][id2]
    return nothing
end


###############################################################################
mutable struct MoleculeFlip <: Action
    i::Int
    j::Int
end

function MonteCarlo.perform_action!(system::Particles, action::MoleculeFlip)
    i, j = action.i, action.j
    spi, spj = system.species[i], system.species[j]
    e₁, e₂ = swap_particle_species!(system, spi, i, spj, j)
    return e₁, e₂
end

function MonteCarlo.invert_action!(action::MoleculeFlip, ::Molecules)
    action.i, action.j = action.j, action.i
    return nothing
end

function reward(::MoleculeFlip, system::Particles)
    return one(typeof(system.temperature))
end

function MonteCarlo.log_proposal_density(action::MoleculeFlip, ::DoubleUniform, parameters, system::Particles)
    return -log(2)
end

function MonteCarlo.sample_action!(action::MoleculeFlip, ::DoubleUniform, parameters, system::Particles, rng)
    moleculei = rand(rng, DiscreteUniform(1, system.N_mol))
    mol_speciesi = system.mol_species[moleculei]
    start_mol, end_mol = mol_speciesi[3], mol_speciesi[4]
    action.i, action.j = sample(rng, start_mol:end_mol, 2; replace=false)
    return nothing
end

###############################################################################

nothing