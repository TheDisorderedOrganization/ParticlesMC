using Arianna
using ParticlesMC
using StaticArrays
using Distributions
using Random
using ComponentArrays

seed = 42
rng = Xoshiro(seed)
# N = 10
# N = 16
N = 44
@assert iseven(N)
NA = N ÷ 2
NB = NA
density = 0.5
d = 2

# temperature = 2.0
# temperature = 1.0

# temperature = 0.5
# temperature = 0.3
# temperature = 0.2
temperature = 0.1 
# temperature = 0.08
# temperature = 0.075
# temperature = 0.07
# temperature = 0.06
# temperature = 0.05
# temperature = 0.04
# temperature = 0.01

# # Analysis run
# M = 1
# steps = 10^6
# burn = 100
# block = [0, 1, 2, 4, 8, 16, 32, 64, 128]
# path = "data/test/particles/SS142D/test/T$temperature/N$N/M$M/steps$steps/seed$seed"

# Dataset run
M = 100
burn = 10^4
# burn = 5 * 10^4
# burn = 10^6
# steps = 100 * burn
# steps = 500 * burn # FIVE TIMES LARGER DATASET
steps = 1000 * burn # TEN TIMES LARGER DATASET
block = [0, burn]
path = "data/datasets/SS142D/T$temperature/N$N/M$M/steps$steps/seed$seed"



box = @SVector fill(typeof(temperature)((N / density)^(1 / d)), d)
position = [[box .* @SVector rand(rng, d) for i in 1:N] for m in 1:M]
species = [shuffle!(rng, vcat(ones(Int, NA), 2ones(Int, NB))) for _ in 1:M]
model = BHHP()
chains = [System(position[m], species[m], density, temperature, model) for m in 1:M]
displacement_policy = SimpleGaussian()
displacement_parameters = ComponentArray(σ=0.065)
pool = (
    Move(Displacement(0, zero(box)), displacement_policy, displacement_parameters, 1.0),
)
callbacks = (callback_energy, callback_acceptance)


sampletimes = build_schedule(steps, burn, block)
algorithm_list = (
    (algorithm=Metropolis, pool=pool, seed=seed, parallel=true, sweepstep=N),
    (algorithm=StoreCallbacks, callbacks=(callback_energy, callback_acceptance), scheduler=sampletimes, store_first=false, store_last=false,),
    (algorithm=StoreTrajectories, scheduler=sampletimes, store_first=false, store_last=false, fmt=XYZ()),
    (algorithm=StoreLastFrames, scheduler=[steps], fmt=XYZ()),
    (algorithm=PrintTimeSteps, scheduler=build_schedule(steps, burn, steps ÷ 10)),
)
simulation = Simulation(chains, algorithm_list, steps; path=path, verbose=true)
run!(simulation)