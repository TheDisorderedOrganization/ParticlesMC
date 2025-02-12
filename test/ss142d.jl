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
model = ParticlesMC.BHHP()
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
sampletimes = scheduler(steps, burn, block)
path = "data/test/particles/SS142D/T$temperature/N$N/M$M/steps$steps/seed$seed"
simulation = Simulation(chains, pools, steps; sweepstep=N, sampletimes=sampletimes, seed=seed, store_trajectory=true, parallel=true, verbose=true, path=path)
callbacks = (callback_energy, callback_acceptance)
## Run the simulation
run!(simulation, callbacks...)