# TEST WITH GERHARD'S CONFIG
using Arianna
using ParticlesMC
using StaticArrays
using Distributions
using Random
using ComponentArrays
using BenchmarkTools
using Profile, PProf
chains_el = load_chains("test/config_0.lmp", args=Dict("temperature" => [0.231], "model" => ["JBB"], "list_type" => "EmptyList"))
chains_ll = load_chains("test/config_0.lmp", args=Dict("temperature" => [0.231], "model" => ["JBB"], "list_type" => "CellList"))

system_el = chains_el[1]
system_ll = chains_ll[1]
# GERHARD: -2.676832
@show mean(system_el.local_energy) / 2
@show mean(system_ll.local_energy) / 2
chains_bkp = deepcopy(chains_el)
chains_ll_bkp = deepcopy(chains_ll)

M = 1
seed = 10
rng = Xoshiro(seed)
sp1 = findall(isequal(1), system_el.species)
sp2 = findall(isequal(2), system_el.species)
sp3 = findall(isequal(3), system_el.species)
NA = length(sp1)
NB = length(sp2)
NC = length(sp3)
steps = 5000
burn = 0
block = [0, 100]

# NO SWAPS
pswap = 0.0
action = Displacement(1, [0.01, 0.01], 0.0)
displacement_policy = SimpleGaussian()
displacement_parameters = ComponentArray(σ=0.05)
pool = (
    Move(Displacement(0, zero(system_ll.box), 0.0), displacement_policy, displacement_parameters, 1 - pswap),
) 
## Empty List
N = length(system_ll)
temperature = system_ll.temperature
chains = deepcopy(chains_bkp)
path = "data/test/particles/KA2D_distribution/N$N/T$temperature/pswap$pswap/M$M"
sampletimes = build_schedule(steps, burn, block)
callbacks = (callback_energy, callback_acceptance)

algorithm_list = (
    (algorithm=Metropolis, pool=pool, seed=seed, parallel=false, sweepstep=length(system_ll)),
    (algorithm=StoreCallbacks, callbacks=(callback_energy, callback_acceptance), scheduler=sampletimes),
    (algorithm=StoreTrajectories, scheduler=sampletimes, fmt=XYZ()),
    (algorithm=StoreLastFrames, scheduler=[steps], fmt=XYZ()),
    (algorithm=PrintTimeSteps, scheduler=build_schedule(steps, burn, steps ÷ 10)),
)
simulation = Simulation(chains_ll_bkp, algorithm_list, steps; path=path, verbose=true)

run!(simulation)



#displacement_action = Displacement(0, zero(box))
#@btime mc_step!(system, displacement_action, displacement_policy, displacement_parameters, rng)

# Profile.clear()
# @profile run!(simulation, callbacks...)
# pprof()



# function test(num, system, action::Action, policy::Policy, parameters::AbstractArray{T}, rng) where {T<:AbstractFloat}
#     for i in 1:num
#         println(i)
#         mc_step!(system, action, policy, parameters, rng)
#     end
# end


# displacement_action = Displacement(0, zero(box))
# test(1, system, displacement_action, displacement_policy, displacement_parameters, rng)
# @code_warntype test(1, system, displacement_action, displacement_policy, displacement_parameters, rng)
# @code_warntype mc_step!(system, displacement_action, displacement_policy, displacement_parameters, rng)

# Profile.clear()
# @profile begin
#     test(1000000, system, displacement_action, displacement_policy, displacement_parameters, rng)
# end
# pprof()

## Linked List
#steps=1000
#sampletimes = scheduler(steps, burn, block)
#chains = deepcopy(chains_ll_bkp)
#path = "data/test/particles/KA2D_distribution_LL/N$N/T$temperature/pswap$pswap/M$M"
#simulation = Simulation(chains, pools, steps; sweepstep=N, sampletimes=sampletimes, seed=seed, parallel=false, verbose=true, path=path)

#Profile.clear()
#@profile run!(simulation, callbacks...)
#pprof()