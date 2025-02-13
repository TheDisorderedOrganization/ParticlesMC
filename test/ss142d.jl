using MonteCarlo
using ParticlesMC
using StaticArrays
using Distributions
using Random
using ComponentArrays

seed = 42
rng = Xoshiro(seed)
NA = 4
NB = 4
N = NA + NB
M = 100
d = 2
temperature = 1.0
density = 0.5
box = @SVector fill(typeof(temperature)((N / density)^(1 / d)), d)
position = [[box .* @SVector rand(rng, d) for i in 1:N] for m in 1:M]
species = [shuffle!(rng, vcat(ones(Int, NA), 2ones(Int, NB))) for _ in 1:M]
model = BHHP()
chains = [System(position[m], species[m], density, temperature, model) for m in 1:M]
displacement_policy = SimpleGaussian()
displacement_parameters = ComponentArray(σ=0.065)
pools = [(
    Move(Displacement(0, zero(box)), displacement_policy, displacement_parameters, 1.0),
) for _ in 1:M]
## Define the simulation struct
steps = 10^6
# burn = 100
# block = [0, 1, 2, 4, 8, 16, 32, 64, 128]
burn = 10^3
block = [0, burn]
callbacks = (callback_energy, callback_acceptance)

path = "data/test/particles/SS142D/T$temperature/N$N/M$M/steps$steps/seed$seed"
sampletimes = build_schedule(steps, burn, block)
algorithm_list = (
    (algorithm=Metropolis, pools=pools, seed=seed, parallel=false, sweepstep=N),
    (algorithm=StoreCallbacks, callbacks=(callback_energy, callback_acceptance), scheduler=sampletimes),
    (algorithm=StoreTrajectories, scheduler=sampletimes),
    (algorithm=StoreLastFrames, scheduler=[steps]),
    (algorithm=PrintTimeSteps, scheduler=build_schedule(steps, burn, steps ÷ 10)),
)
simulation = Simulation(chains, algorithm_list, steps; path=path, verbose=true)
run!(simulation)