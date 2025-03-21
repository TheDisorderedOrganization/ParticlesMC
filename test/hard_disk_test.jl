using Arianna
using Arianna.PolicyGuided
using ParticlesMC
using StaticArrays
using Distributions
using Random
using ComponentArrays

seed = 42
rng = Xoshiro(seed)
NA = 108
NB = 108
N = NA + NB
M = 1
d = 2
temperature = 1.0
density = 0.2
pressure = 10.0
box = @SVector fill(typeof(temperature)((N / density)^(1 / d)), d)
position = [[box .* @SVector rand(rng, d) for i in 1:N] for m in 1:M]
species = [shuffle!(rng, vcat(ones(NA), 1.4 * ones(NB))) for _ in 1:M]

chains = [System(position[m], species[m], density, temperature, HardCore(); list_type=LinkedList) for m in 1:M]
displacement_policy = SimpleGaussian()
displacement_parameters = ComponentArray(σ=0.5)
barostat_policy = SimpleGaussian()
barostat_parameters = ComponentArray(σ=0.5)
# pool = (
#     Move(Displacement(0, zero(box)), displacement_policy, displacement_parameters, 1.0),
# )
pool = (
    Move(Displacement(0, zero(box)), displacement_policy, displacement_parameters, 1 - 1 / N),
    Move(Barostat(pressure, 0.0), barostat_policy, barostat_parameters, 1 / N),
)
optimisers = (BLAPG(1e-5, 1e-4), Static())
steps = 10^5
burn = 0
block = [0, 1, 2, 4, 8, 16, 32, 64, 128]
# burn = 10^3
# block = [0, burn]

pgmc_start = 10^4
estimator_scheduler = build_schedule(steps, pgmc_start, 1)
learner_scheduler = build_schedule(steps, pgmc_start, 10)

path = "data/test/particles/HardDisks/"
sampletimes = build_schedule(steps, burn, block)
algorithm_list = (
    (algorithm=Metropolis, pool=pool, seed=seed, parallel=false, sweepstep=N),
    (algorithm=PolicyGradientEstimator, dependencies=(Metropolis,), optimisers=optimisers, q_batch_size=20, parallel=false, scheduler=estimator_scheduler),
    (algorithm=PolicyGradientUpdate, dependencies=(PolicyGradientEstimator,), scheduler=learner_scheduler),
    (algorithm=StoreCallbacks, callbacks=(callback_acceptance,), scheduler=sampletimes),
    (algorithm=StoreTrajectories, scheduler=sampletimes, fmt=XYZ()),
    (algorithm=StoreLastFrames, scheduler=[steps], fmt=XYZ()),
    (algorithm=StoreParameters, dependencies=(Metropolis,), scheduler=sampletimes),
    (algorithm=PrintTimeSteps, scheduler=build_schedule(steps, burn, steps ÷ 10)),
)
simulation = Simulation(chains, algorithm_list, steps; path=path, verbose=true)
run!(simulation)