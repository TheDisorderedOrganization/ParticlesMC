"""
ParticlesMC: Monte Carlo simulation framework for particle systems.

Provides core types, utilities, and the `particlesmc` command implemented with Comonicon.
Exports commonly-used types (e.g., `Particles`, `Model`) and helper functions for simulation control, I/O, and moves.
"""
module ParticlesMC

using Arianna, StaticArrays
using Comonicon, TOML
using Comonicon: @main

export Particles
abstract type Particles <: AriannaSystem end


include("utils.jl")
include("neighbours.jl")
include("models.jl")
include("molecules.jl")
include("atoms.jl")
include("moves.jl")

"""Return the position of particle `i` in `system`.

# Arguments
- `system::Particles`: the particle system
- `i::Int`: particle index

# Returns
- Coordinates of particle `i` (e.g., an `SVector` or array).
"""
get_position(system::Particles, i::Int) = @inbounds system.position[i]

"""Return the species index (type) of particle `i`.

# Arguments
- `system::Particles`: the particle system
- `i::Int`: particle index

# Returns
- `Int`: species identifier of particle `i`.
"""
get_species(system::Particles, i::Int) = @inbounds system.species[i]

"""Return the interaction `Model` between the species of particles `i` and `j`.

# Arguments
- `system::Particles`: the particle system
- `i::Int`, `j::Int`: particle indices

# Returns
- `Model` object or callable describing pair interactions for the two species.
"""
get_model(system::Particles, i::Int, j::Int) = @inbounds system.model_matrix[get_species(system, i), get_species(system, j)]

"""Return the simulation box of `system`.

# Returns
- Box description (usually vector or struct) representing periodic box extents.
"""
get_box(system::Particles) = system.box

"""Return the neighbour list of `system`.

# Returns
- The neighbour list object used for pair evaluations (e.g., `NeighbourList`, `LinkedList`).
"""
get_neighbour_list(system::Particles) = system.neighbour_list

"""Return the number of particles in `system`.

Overloads `Base.length` for `Particles`.
"""
Base.length(system::Particles) = system.N

"""Return a proper index range for `system`.

Overloads `Base.eachindex` to allow fast indexing.
"""
Base.eachindex(system::Particles) = Base.OneTo(length(system))

"""Return the (position, species) tuple for atom `i`.

This overload supports indexing `atoms[i]` to get coordinates and species.
"""
Base.getindex(system::Atoms, i::Int) = system.position[i], system.species[i]

"""Iterate over `Atoms` or `Molecules` returning the position and next state.

Conforms to Julia iterator interface; yields the position of the current index.
"""
function Base.iterate(system::Union{Atoms, Molecules}, state=1)
    state > length(system) && return nothing  # Stop iteration
    return (system.position[state], state + 1)  # Return element & next state
end

"""Compute the energy contribution of particle `i` in `system` using its neighbour list.

If a neighbour list is provided it will be used for the evaluation.
"""
function compute_energy_particle(system::Particles, i::Int)
    return compute_energy_particle(system, i, system.neighbour_list)
end

"""Compute the energy contribution for each particle index in `ids`.

Returns an array with per-particle energies by mapping `compute_energy_particle`.
"""
function compute_energy_particle(system::Particles, ids::AbstractVector)
    return map(i -> compute_energy_particle(system, i), ids)
end


export callback_energy
#export nearest_image_distance
export Model, GeneralKG, JBB, BHHP, SoftSpheres, KobAndersen, Trimer
export NeighbourList, LinkedList, CellList, EmptyList
export Atoms, Molecules
export Displacement, DiscreteSwap, MoleculeFlip
export fold_back, System
export SimpleGaussian, DoubleUniform, EnergyBias
export sample_action!, log_proposal_density, reward, invert_action!, delta_log_target_density
export perform_action!, revert_action!
include("IO/IO.jl")
using .IO: XYZ, EXYZ, LAMMPS, load_configuration, load_chains
export XYZ, EXYZ, LAMMPS, load_configuration, load_chains


"""
ParticlesMC implemented in Comonicon.

# Arguments

- `params`: Path to the TOML parameter file.
"""
@main function particlesmc(params::String)
    if !isfile(params)
        error("Parameter file '$params' does not exist in the current path.")
    end
    params = TOML.parsefile(params)

    # Extract system parameters
    system = params["system"]
    temperature = system["temperature"]
    density = system["density"]
    config = system["config"]
    model = get(system, "model", nothing)
    if model === nothing
        model = params["model"]
    end  # optional field
    if !isfile(config)
        error("Configuration file '$config' does not exist in the current path.")
    end
    list_type = get(system, "list_type", "LinkedList")  # optional field
    bonds = get(system, "bonds", nothing)  
    
    # Extract simulation parameters
    sim = params["simulation"]
    steps = sim["steps"]
    burn = get(sim, "burn", 0)
    seed = sim["seed"]
    parallel = sim["parallel"]
    output_path = get(sim, "output_path", "./")

    # Setup RNG and basic variables

    # optional field

    if bonds !== nothing
        chains = load_chains(config, args=Dict(
            "temperature" => [temperature],
            "density" => [density],
            "model" => [model],
            "list_type" => list_type,
            "bonds" => bonds,
        ))
    else    
        chains = load_chains(config, args=Dict(
            "temperature" => [temperature],
            "density" => [density],
            "model" => [model],
            "list_type" => list_type,
        ))
    end
    algorithm_list = []
    # Setup moves
    pool = []
    for move in sim["move"]
        prob = move["probability"]
        policy = move["policy"]
        action = move["action"]
        parameters = get(move, "parameters", Dict())
        param_obj = ComponentArray()

        # Create action object
        if action == "Displacement"
            action_obj = Displacement(0, zero(chains[1].box), 0.0)
            if "sigma" in keys(parameters)
                param_obj = ComponentArray(σ = parameters["sigma"])
            else
                error("Missing parameter 'sigma' for action: $action")
            end
            if policy == "SimpleGaussian"
                policy_obj = SimpleGaussian()
            else
                error("Unsupported policy: $policy for action: $action")
            end
        elseif action == "MoleculeFlip"
            action_obj = MoleculeFlip(0, 0, 0.0)
            param_obj = Vector{Float64}()
            if policy == "DoubleUniform"
                policy_obj = DoubleUniform()
            else
                error("Unsupported policy: $policy for action: $action")
            end
        else
            error("Unsupported action: $action") 
        end
        # Build move
        move_obj = Move(action_obj, policy_obj, param_obj, prob)
        push!(pool, move_obj)
    end
    push!(algorithm_list, (algorithm=Metropolis, pool=pool, seed=seed, parallel=parallel, sweepstep=length(chains[1])))

    # Setup outputs
    for output in sim["output"]
        alg = output["algorithm"]
        scheduler_params = output["scheduler_params"]
        callbacks = get(output, "callbacks", [])
        fmt = get(output, "fmt", "XYZ")
        interval = scheduler_params["linear_interval"]
        if "log_base" in keys(scheduler_params)
            block = build_schedule(interval, 0,  2.0)
            sched = build_schedule(steps, burn, block)
        else
            sched = build_schedule(steps, burn, interval)
        end
        if alg == "StoreCallbacks"
            callbacks = map(c -> eval(Meta.parse("callback_$c")), callbacks)
            algorithm = (
                algorithm = eval(Meta.parse(alg)),
                callbacks = callbacks,
                scheduler = sched,
            )
        elseif alg == "StoreTrajectories" || alg == "StoreLastFrames"
             algorithm = (
                algorithm = eval(Meta.parse(alg)),
                scheduler = sched,
                fmt = eval(Meta.parse("$(fmt)()")),
            )
        elseif alg == "PrintTimeSteps"
            algorithm = (
                algorithm = eval(Meta.parse(alg)),
                scheduler = build_schedule(steps, burn, steps ÷ 10),
            )
        else
            error("Unsupported output algorithm: $alg")
        end
        push!(algorithm_list, algorithm)    
    end
    M=1
    path = joinpath(output_path)
    simulation = Simulation(chains, algorithm_list, steps; path=path, verbose=true)
    
    # Run the simulation
    run!(simulation)

end 

end