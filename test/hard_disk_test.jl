using Arianna
using ParticlesMC
using StaticArrays
using Distributions
using Random
using ComponentArrays

seed = 42
rng = Xoshiro(seed)
NA = 5
NB = 5
N = NA + NB
M = 1
d = 2
temperature = 1.0
density = 0.2
pressure = 1.0
box = @SVector fill(typeof(temperature)((N / density)^(1 / d)), d)
position = [[box .* @SVector rand(rng, d) for i in 1:N] for m in 1:M]
species = [shuffle!(rng, vcat(ones(NA), 1.4 * ones(NB))) for _ in 1:M]

#chains = [HardSpheres(position[m], species[m], density, temperature; list_type=EmptyList) for m in 1:M]
chains = [HardSpheres(position[m], species[m], density, temperature; list_type=LinkedList) for m in 1:M]
displacement_policy = SimpleGaussian()
displacement_parameters = ComponentArray(σ=0.5)
pool = (
    Move(Displacement(0, zero(box)), displacement_policy, displacement_parameters, 1.0),
)
# pool = (
#     Move(Displacement(0, zero(box)), displacement_policy, displacement_parameters, 1 - 1 / N),
#     Move(Barostat(pressure, 0.0), displacement_policy, displacement_parameters, 1 / N),
# )
steps = 10^4
burn = 0
block = [0, 1, 2, 4, 8, 16, 32, 64, 128]
# burn = 10^3
# block = [0, burn]

path = "data/test/particles/HardDisks/"
sampletimes = build_schedule(steps, burn, block)
algorithm_list = (
    (algorithm=Metropolis, pool=pool, seed=seed, parallel=false, sweepstep=N),
    (algorithm=StoreCallbacks, callbacks=(callback_acceptance,), scheduler=sampletimes),
    (algorithm=StoreTrajectories, scheduler=sampletimes, fmt=XYZ()),
    (algorithm=StoreLastFrames, scheduler=[steps], fmt=XYZ()),
    (algorithm=PrintTimeSteps, scheduler=build_schedule(steps, burn, steps ÷ 10)),
)
simulation = Simulation(chains, algorithm_list, steps; path=path, verbose=true)
run!(simulation)