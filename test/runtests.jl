using Arianna
using ParticlesMC
using Test
using StaticArrays
using Distributions
using ComponentArrays
using DelimitedFiles
using Pkg

@testset "Running from CLI" begin
    Pkg.build("ParticlesMC")
    path_sep = Sys.iswindows() ? ";" : ":"
    julia_bin = expanduser("~/.julia/bin")
    ENV["PATH"] = ENV["PATH"] * path_sep * julia_bin
    @test success(`bash -c "command -v particlesmc"`)
    @test success(`particlesmc params.toml`)
end 


@testset "Potential energy test" begin
    # Test inital configuration
    chains_el = load_chains("config_0.exyz", args=Dict("temperature" => [0.231], "model" => ["JBB"], "list_type" => "EmptyList"))
    chains_ll = load_chains("config_0.lmp", args=Dict("temperature" => [0.231], "model" => ["JBB"], "list_type" => "LinkedList"))
    system_el = chains_el[1]
    system_ll = chains_ll[1]
    @test system_el.N == system_ll.N
    @test system_el.d == system_ll.d
    @test system_el.temperature == system_ll.temperature
    @test all.(system_el.position == system_ll.position)
    @test all.(system_el.species == system_ll.species)
    energy_el = system_el.energy[1] / length(system_el)
    energy_ll = system_ll.energy[1] / length(system_ll)
    @test isapprox(energy_el, -2.676832, atol=1e-6)
    @test isapprox(energy_ll, -2.676832, atol=1e-6)

    # Test simulation energy
    M = 1
    seed = 10
    sp1, sp2, sp3 = findall(isequal(1), system_el.species), findall(isequal(2), system_el.species), findall(isequal(3), system_el.species)
    NA, NB, NC = length(sp1), length(sp2), length(sp3)
    steps = 100
    burn = 0
    block = [0, 1, 2, 4, 8]
    sampletimes = build_schedule(steps, burn, block)
    callbacks = (callback_energy, callback_acceptance)

    # NO SWAPS
    pswap = 0.0
    displacement_policy = SimpleGaussian()
    displacement_parameters = ComponentArray(σ=0.05)
    pool = (
        Move(Displacement(0, zero(system_el.box), 0.0), displacement_policy, displacement_parameters, 1 - pswap),
    )
    algorithm_list = (
    (algorithm=Metropolis, pool=pool, seed=seed, parallel=false, sweepstep=system_el.N),
    (algorithm=StoreCallbacks, callbacks=(callback_energy, callback_acceptance), scheduler=sampletimes),
    (algorithm=StoreTrajectories, scheduler=sampletimes, fmt=EXYZ()),
    (algorithm=StoreLastFrames, scheduler=[steps], fmt=LAMMPS()),
    (algorithm=PrintTimeSteps, scheduler=build_schedule(steps, burn, steps ÷ 10)),
    )
    ## Empty List simulation
    chains_el = [deepcopy(system_el)]
    path_el = "data/noswap/empty_list/"
    simulation = Simulation(chains_el, algorithm_list, steps; path=path_el, verbose=true)
    run!(simulation)
    
    ## Linked List simulation
    chains_ll = [deepcopy(system_ll)]
    path_ll = "data/noswap/linked_list/"
    simulation = Simulation(chains_ll, algorithm_list, steps; path=path_ll, verbose=true)
    run!(simulation)

    ## Read energy data and compare
    path_energy_el = joinpath(path_el, "energy.dat")
    path_energy_ll = joinpath(path_ll, "energy.dat")
    energy_el= readdlm(path_energy_el)[:, 2]
    energy_ll = readdlm(path_energy_ll)[:, 2]
    @test isapprox(energy_el, energy_ll, atol=1e-6)

    # SWAPS
    pswap = 0.8
    displacement_policy = SimpleGaussian()
    displacement_parameters = ComponentArray(σ=0.05)
    swap_policy = DoubleUniform()
    swap_parameters = Vector{Float64}()
    pool = (
    Move(Displacement(0, zero(system_el.box), 0.0), displacement_policy, displacement_parameters, 1 - pswap),
    Move(DiscreteSwap(0, 0, (1, 3), (NA, NC), 0.0), swap_policy, swap_parameters, pswap / 2),
    Move(DiscreteSwap(0, 0, (2, 3), (NB, NC), 0.0), swap_policy, swap_parameters, pswap / 2),
    )
    algorithm_list = (
        (algorithm=Metropolis, pool=pool, seed=seed, parallel=false, sweepstep=system_el.N),
        (algorithm=StoreCallbacks, callbacks=(callback_energy, callback_acceptance), scheduler=sampletimes, fmt=XYZ()),
        (algorithm=StoreTrajectories, scheduler=sampletimes, fmt=XYZ()),
        (algorithm=StoreLastFrames, scheduler=[steps], fmt=XYZ()),
        (algorithm=PrintTimeSteps, scheduler=build_schedule(steps, burn, steps ÷ 10), fmt=XYZ()),
        )
    ## Empty List simulation
    chains_el = [deepcopy(system_el)]
    path_el = "data/swap/empty_list/"
    simulation = Simulation(chains_el, algorithm_list, steps; path=path_el, verbose=true)
    run!(simulation)
    
    ## Linked List simulation
    chains_ll = [deepcopy(system_ll)]
    path_ll = "data/swap/linked_list/"
    simulation = Simulation(chains_ll, algorithm_list, steps; path=path_ll, verbose=true)
    run!(simulation)

    ## Read energy data and compare
    path_energy_el = joinpath(path_el, "energy.dat")
    path_energy_ll = joinpath(path_ll, "energy.dat")
    energy_el= readdlm(path_energy_el)[:, 2]
    energy_ll = readdlm(path_energy_ll)[:, 2]
    @test isapprox(energy_el, energy_ll, atol=1e-6)


end

@testset "Molecule potential energy test" begin
    # Test inital configuration
    chains_el = load_chains("molecule.exyz", args=Dict("temperature" => [2.0], "model" => ["Trimer"], "list_type" => "EmptyList"))
    chains_ll = load_chains("molecule.xyz", args=Dict("temperature" => [2.0], "model" => ["Trimer"], "list_type" => "LinkedList"))
    system_el = chains_el[1]
    system_ll = chains_ll[1]
    @test system_el.N == system_ll.N
    @test system_el.d == system_ll.d
    @test system_el.temperature == system_ll.temperature
    @test all.(system_el.position == system_ll.position)
    @test all.(system_el.species == system_ll.species)
    @test all.(system_el.bonds == system_ll.bonds)
    energy_el = system_el.energy[1] / length(system_el)
    energy_ll = system_ll.energy[1] / length(system_ll)
    @test isapprox(energy_el, 25.65865662277199, atol=1e-6)
    @test isapprox(energy_ll, 25.65865662277199, atol=1e-6)

    # Test simulation energy
    M = 1
    seed = 10
    steps = 100
    burn = 0
    block = [0, 1, 2, 4, 8]
    sampletimes = build_schedule(steps, burn, block)
    callbacks = (callback_energy, callback_acceptance)

    # NO SWAPS
    pswap = 0.0
    displacement_policy = SimpleGaussian()
    displacement_parameters = ComponentArray(σ=0.05)
    pool = (
        Move(Displacement(0, zero(system_el.box), 0.0), displacement_policy, displacement_parameters, 1 - pswap),
    )
    algorithm_list = (
    (algorithm=Metropolis, pool=pool, seed=seed, parallel=false, sweepstep=system_el.N),
    (algorithm=StoreCallbacks, callbacks=(callback_energy, callback_acceptance), scheduler=sampletimes),
    (algorithm=StoreTrajectories, scheduler=sampletimes, fmt=EXYZ()),
    (algorithm=StoreLastFrames, scheduler=[steps], fmt=EXYZ()),
    (algorithm=PrintTimeSteps, scheduler=build_schedule(steps, burn, steps ÷ 10)),
    )
    ## Empty List simulation
    chains_el = [deepcopy(system_el)]
    path_el = "data/noswap/empty_list/"
    simulation = Simulation(chains_el, algorithm_list, steps; path=path_el, verbose=true)
    run!(simulation)
    
    ## Linked List simulation
    chains_ll = [deepcopy(system_ll)]
    path_ll = "data/noswap/linked_list/"
    simulation = Simulation(chains_ll, algorithm_list, steps; path=path_ll, verbose=true)
    run!(simulation)

    ## Read energy data and compare
    path_energy_el = joinpath(path_el, "energy.dat")
    path_energy_ll = joinpath(path_ll, "energy.dat")
    energy_el= readdlm(path_energy_el)[:, 2]
    energy_ll = readdlm(path_energy_ll)[:, 2]
    @test isapprox(energy_el, energy_ll, atol=1e-6)

end


