using MonteCarlo
using MonteCarlo.PolicyGuided
using ParticlesMC
using StaticArrays
using Distributions
using Random
using ComponentArrays

seed = 42
rng = Xoshiro(seed)
NA = 20
NB = 11
NC = 12
N = NA + NB + NC
M = 10
d = 2
temperature = 0.5
density = 1.1920748468939728
box = @SVector fill(typeof(temperature)((N / density)^(1 / d)), d)
position = [[box .* @SVector rand(rng, d) for i in 1:N] for m in 1:M]
species = [shuffle!(rng, vcat(ones(Int, NA), 2ones(Int, NB), 3ones(Int, NC))) for _ in 1:M]
model = JBB()
chains = [System(position[m], species[m], density, temperature, model) for m in 1:M]
pswap = 0.2
displacement_policy = SimpleGaussian()
displacement_parameters = ComponentArray(σ=0.05)
swap_AC_policy = EnergyBias()
swap_BC_policy = EnergyBias()
swap_AC_parameters = ComponentArray(θ₁=0.0, θ₂=0.0)
swap_BC_parameters = ComponentArray(θ₁=0.0, θ₂=0.0)
pool = (
    Move(Displacement(0, zero(box)), displacement_policy, displacement_parameters, 1 - pswap),
    Move(DiscreteSwap(0, 0, (1, 3), (NA, NC)), swap_AC_policy, swap_AC_parameters, pswap / 2),
    Move(DiscreteSwap(0, 0, (2, 3), (NB, NC)), swap_BC_policy, swap_BC_parameters, pswap / 2),
)
optimisers = (VPG(1e-3), BLANPG(1e-6, 1e-6), BLANPG(1e-6, 1e-6))
# optimisers = (VPG(1e-5), VPG(1e-2), VPG(1e-2))
steps = 50000
burn = 10000
block = [0, 1, 2, 4, 8, 16, 32, 64, 128]
sampletimes = build_schedule(steps, burn, block)
callbacks = (callback_energy, callback_acceptance)

path = "data/test/pgmc/KA2D/T$temperature/N$N/M$M/steps$steps/seed$seed"
sampletimes = build_schedule(steps, burn, block)
algorithm_list = (
    (algorithm=Metropolis, pool=pool, seed=seed, parallel=true, sweepstep=N),
    (algorithm=PolicyGradientEstimator, dependencies=(Metropolis,), optimisers=optimisers, q_batch_size=10, parallel=true),
    (algorithm=PolicyGradientUpdate, dependencies=(PolicyGradientEstimator,), scheduler=build_schedule(steps, burn, 2)),
    (algorithm=StoreCallbacks, callbacks=(callback_energy, callback_acceptance), scheduler=sampletimes),
    (algorithm=StoreTrajectories, scheduler=sampletimes),
    (algorithm=StoreParameters, dependencies=(Metropolis,), scheduler=sampletimes),
    (algorithm=StoreLastFrames, scheduler=[steps]),
    (algorithm=PrintTimeSteps, scheduler=build_schedule(steps, burn, steps ÷ 10)),
)
simulation = Simulation(chains, algorithm_list, steps; path=path, verbose=true)
run!(simulation)

