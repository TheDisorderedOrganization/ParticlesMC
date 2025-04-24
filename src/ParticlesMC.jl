module ParticlesMC

using Arianna, StaticArrays

export Particles
abstract type Particles <: AriannaSystem end


include("utils.jl")
include("neighbours.jl")
include("models.jl")
include("molecules.jl")
include("atoms.jl")
include("moves.jl")

get_position(system::Particles, i::Int) = @inbounds system.position[i]
get_species(system::Particles, i::Int) = @inbounds system.species[i]
get_model(system::Particles, i::Int, j::Int) = @inbounds system.model_matrix[get_species(system, i), get_species(system, j)]
get_local_energy(system::Particles, i::Int) = @inbounds system.local_energy[i]
get_box(system::Particles) = system.box
get_neighbour_list(system::Particles) = system.neighbour_list
Base.length(system::Particles) = system.N
Base.eachindex(system::Particles) = Base.OneTo(length(system))
Base.getindex(system::Atoms, i::Int) = system.position[i], system.species[i], system.local_energy[i]

function Base.iterate(system::Union{Atoms, Molecules}, state=1)
    state > length(system) && return nothing  # Stop iteration
    return (system.position[state], state + 1)  # Return element & next state
end

function compute_energy_particle(system::Particles, i)
    return compute_energy_particle(system, i, system.neighbour_list)
end


export callback_energy
export nearest_image_distance
export Model
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

end