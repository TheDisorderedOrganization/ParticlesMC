using Arianna
using ParticlesMC
using Distributions
using Random
using StaticArrays
using ComponentArrays

seed = 42
rng = Xoshiro(seed)
N = 3000
M = 1
Length = 3
d = 3
temperature = 2.0
density = 1.2

function create_bond_matrix(N::Int)
    # Create a vector to store the SVectors, each containing a pair of integers
    matrix = Vector{Vector{Int}}()
    # Populate the matrix with pairs according to the specified pattern
    for i in 2:3:N
        push!(matrix, [i, i + 1])
        push!(matrix, [i - 1, i + 1])
        push!(matrix, [i - 1, i])
    end
    return matrix
end

chains = load_chains("test/molecule.xyz", args=Dict("temperature" => [temperature], "model" => ["Trimer"], "list_type" => "LinkedList", "bonds" => create_bond_matrix(N)))
pswap = 0.2
displacement_policy = SimpleGaussian()
displacement_parameters = ComponentArray(σ=0.05)
pool = (
    Move(Displacement(0, zero(chains[1].box), 0.0), displacement_policy, displacement_parameters, 1.0),
)
## Define the simulation struct
steps = 1000
burn = 0
block = [0, 1000]
sampletimes = build_schedule(steps, burn, block)
callbacks = (callback_energy, callback_acceptance)

path = "data/test/particles/Molecules/T$temperature/N$N/M$M/seed$seed"
algorithm_list = (
    (algorithm=Metropolis, pool=pool, seed=seed, parallel=false, sweepstep=N),
    (algorithm=StoreCallbacks, callbacks=(energy,), scheduler=sampletimes),
    (algorithm=StoreAcceptance, dependencies=(Metropolis,), scheduler=sampletimes),    (algorithm=StoreTrajectories, scheduler=sampletimes, fmt=XYZ()),
    (algorithm=StoreLastFrames, scheduler=[steps], fmt=XYZ()),
    (algorithm=PrintTimeSteps, scheduler=build_schedule(steps, burn, steps ÷ 10), fmt=XYZ())
    )
simulation = Simulation(chains, algorithm_list, steps; path=path, verbose=true)
run!(simulation)