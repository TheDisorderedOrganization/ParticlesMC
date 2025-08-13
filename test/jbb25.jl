using Arianna
using ParticlesMC
using StaticArrays
using Distributions
using Random
using ComponentArrays

seed = 42
rng = Xoshiro(seed)
# NA = 5
# NB = 3
# NC = 3
NA = 20
NB = 12
NC = 12
N = NA + NB + NC
d = 2


# temperature = 10.0 # WTF
# temperature = 5.0 # very very high
# temperature = 2.0 # very high
# temperature = 1.0 # high
# temperature = 0.5 # base
# temperature = 0.32 # Gerhard's highest
temperature = 0.2 # Gerhard's lowest NF

# pswap = 0.0
pswap = 0.2


# # TEST
# M = 1
# steps = 2 * 10^6
# burn = 10^6
# block = [0, 1, 2, 4, 8, 16, 32, 64, 128]
# path = "data/test/particles/JBB25/dynamics/pswap$pswap/T$temperature/N$N/M$M/steps$steps/seed$seed"

# DATASET
M = 100
# burn = 10^4
burn = 5 * 10^4
steps = 100 * burn
block = [0, burn]
path = "data/datasets/JBB25/T$temperature/N$N/M$M/steps$steps/seed$seed"

sampletimes = build_schedule(steps, burn, block)
density = 1.1920748468939728
box = @SVector fill(typeof(temperature)((N / density)^(1 / d)), d)
position = [[box .* @SVector rand(rng, d) for i in 1:N] for m in 1:M]
species = [shuffle!(rng, vcat(ones(Int, NA), 2ones(Int, NB), 3ones(Int, NC))) for _ in 1:M]
model = JBB()
chains = [System(position[m], species[m], density, temperature, model) for m in 1:M]
displacement_policy = SimpleGaussian()
displacement_parameters = ComponentArray(σ=0.05)
swap_policy = DoubleUniform()
swap_parameters = Vector{Float64}()

if iszero(pswap)
    pool = (
        Move(Displacement(0, zero(chains[1].box)), displacement_policy, displacement_parameters, 1.0),
    )
else
    pool = (
        Move(Displacement(0, zero(chains[1].box)), displacement_policy, displacement_parameters, 1 - pswap),
        Move(DiscreteSwap(0, 0, (1, 3), (NA, NC)), swap_policy, swap_parameters, pswap / 2),
        Move(DiscreteSwap(0, 0, (2, 3), (NB, NC)), swap_policy, swap_parameters, pswap / 2),
    )
end

callbacks = (callback_energy, callback_acceptance)
algorithm_list = (
    (algorithm=Metropolis, pool=pool, seed=seed, parallel=true, sweepstep=N),
    (algorithm=StoreCallbacks, callbacks=(callback_energy, callback_acceptance), scheduler=sampletimes, store_first=false, store_last=false,),
    (algorithm=StoreTrajectories, scheduler=sampletimes, store_first=false, store_last=false, fmt=XYZ()),
    (algorithm=StoreLastFrames, scheduler=[steps], fmt=XYZ()),
    (algorithm=PrintTimeSteps, scheduler=build_schedule(steps, burn, steps ÷ 10)),
)
simulation = Simulation(chains, algorithm_list, steps; path=path, verbose=true)
run!(simulation)