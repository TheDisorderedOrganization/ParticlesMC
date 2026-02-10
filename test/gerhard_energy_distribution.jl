# TEST WITH GERHARD'S CONFIG
using Arianna
using ParticlesMC
using StaticArrays
using Distributions
using Random
using ComponentArrays
using Profile
chains = load_chains("test/config_0.xyz", args=Dict("temperature" => [0.231], "model" => ["JBB"], "list_type" => "EmptyList"))
chains_ll = load_chains("test/config_0.lmp", args=Dict("temperature" => [0.231], "model" => ["JBB"], "list_type" => "LinkedList"))
system = chains[1]
system_ll = chains_ll[1]
# GERHARD: -2.676832
@show system.energy[1] / length(system)
@show system_ll.energy[1] / length(system_ll)
chains_bkp = deepcopy(chains)
chains_ll_bkp = deepcopy(chains_ll)

M = 1
seed = 10
rng = Xoshiro(seed)
sp1 = findall(isequal(1), system.species)
sp2 = findall(isequal(2), system.species)
sp3 = findall(isequal(3), system.species)
NA = length(sp1)
NB = length(sp2)
NC = length(sp3)
steps = 1000
burn = 0
block = [0, 10]
sampletimes = build_schedule(steps, burn, block)

# NO SWAPS
pswap = 0.0
displacement_policy = SimpleGaussian()
displacement_parameters = ComponentArray(σ=0.05)
pool = (
    Move(Displacement(0, zero(system.box), 0.0), displacement_policy, displacement_parameters, 1 - pswap),
)

algorithm_list = (
    (algorithm=Metropolis, pool=pool, seed=seed, parallel=false, sweepstep=system.N),
    (algorithm=StoreCallbacks, callbacks=(energy,), scheduler=sampletimes),
    (algorithm=StoreAcceptance, dependencies=(Metropolis,), scheduler=sampletimes),
    (algorithm=StoreTrajectories, scheduler=sampletimes, fmt=EXYZ()),
    (algorithm=StoreLastFrames, scheduler=[steps], fmt=EXYZ()),
    (algorithm=PrintTimeSteps, scheduler=build_schedule(steps, burn, steps ÷ 10)),
)

## Empty List
chains = deepcopy(chains_bkp)
path = "data/test/particles/KA2D_distribution/N$(system.N)/T$(system.temperature)/pswap$pswap/M$M"
simulation = Simulation(chains, algorithm_list, steps; path=path, verbose=true)
run!(simulation)

## Linked List
chains = deepcopy(chains_ll_bkp)
path = "data/test/particles/KA2D_distribution_LL/N$(system.N)/T$(system.temperature)/pswap$pswap/M$M"
simulation = Simulation(chains, algorithm_list, steps; path=path, verbose=true)
run!(simulation)

# SWAPS
pswap = 0.2
displacement_policy = SimpleGaussian()
displacement_parameters = ComponentArray(σ=0.05)
swap_policy = DoubleUniform()
swap_parameters = Vector{Float64}()
pool = (
    Move(Displacement(0, zero(system.box), 0.0), displacement_policy, displacement_parameters, 1 - pswap),
    Move(DiscreteSwap(0, 0, (1, 3), (NA, NC), 0.0), swap_policy, swap_parameters, pswap / 2),
    Move(DiscreteSwap(0, 0, (2, 3), (NB, NC), 0.0), swap_policy, swap_parameters, pswap / 2),
)
algorithm_list = (
    (algorithm=Metropolis, pool=pool, seed=seed, parallel=false, sweepstep=system.N),
    (algorithm=StoreCallbacks, callbacks=(energy,), scheduler=sampletimes),
    (algorithm=StoreAcceptance, dependencies=(Metropolis,), scheduler=sampletimes),
    (algorithm=StoreTrajectories, scheduler=sampletimes, fmt=LAMMPS()),
    (algorithm=StoreLastFrames, scheduler=[steps], fmt=EXYZ()),
    (algorithm=PrintTimeSteps, scheduler=build_schedule(steps, burn, steps ÷ 10)),
)
## Empty List
chains = deepcopy(chains_bkp)
path = "data/test/particles/KA2D_distribution/N$(system.N)/T$(system.temperature)/pswap$pswap/M$M"
simulation = Simulation(chains, algorithm_list, steps; path=path, verbose=true)
run!(simulation)

## Linked List
chains = deepcopy(chains_ll_bkp)
path = "data/test/particles/KA2D_distribution_LL/N$(system.N)/T$(system.temperature)/pswap$pswap/M$M"
simulation = Simulation(chains, algorithm_list, steps; path=path, verbose=true)
run!(simulation)


# EB SWAPS
pswap = 0.1
displacement_policy = SimpleGaussian()
displacement_parameters = ComponentArray(σ=0.05)
swap_AC_policy = EnergyBias()
swap_BC_policy = EnergyBias()
swap_AC_parameters = ComponentArray(θ₁=1.0, θ₂=0.5)
swap_BC_parameters = ComponentArray(θ₁=0.5, θ₂=4.0)
pool = (
    Move(Displacement(0, zero(system.box), 0.0), displacement_policy, displacement_parameters, 1 - pswap),
    Move(DiscreteSwap(0, 0, (1, 3), (NA, NC), 0.0), swap_AC_policy, swap_AC_parameters, pswap / 2),
    Move(DiscreteSwap(0, 0, (2, 3), (NB, NC), 0.0), swap_BC_policy, swap_BC_parameters, pswap / 2),
)
algorithm_list = (
    (algorithm=Metropolis, pool=pool, seed=seed, parallel=false, sweepstep=system.N),
    (algorithm=StoreCallbacks, callbacks=(energy,), scheduler=sampletimes),
    (algorithm=StoreAcceptance, dependencies=(Metropolis,), scheduler=sampletimes),
    (algorithm=StoreTrajectories, scheduler=sampletimes, fmt=LAMMPS()),
    (algorithm=StoreLastFrames, scheduler=[steps], fmt=LAMMPS()),
    (algorithm=PrintTimeSteps, scheduler=build_schedule(steps, burn, steps ÷ 10)),
)
## Empty List
chains = deepcopy(chains_bkp)
path = "data/test/particles/KA2D_distribution/N$(system.N)/T$(system.temperature)/ebswap$pswap/M$M"
simulation = Simulation(chains, algorithm_list, steps; path=path, verbose=true)
run!(simulation)

## Linked List
chains = deepcopy(chains_ll_bkp)
path = "data/test/particles/KA2D_distribution_LL/N$(system.N)/T$(system.temperature)/ebswap$pswap/M$M"
simulation = Simulation(chains, algorithm_list, steps; path=path, verbose=true)
run!(simulation)