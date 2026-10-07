using Arianna
using ParticlesMC
using Test
using StaticArrays
using Distributions
using ComponentArrays
using DelimitedFiles
using LinearAlgebra
using Pkg

@testset "Molecular orientation" begin
    box = SVector(10.0, 10.0, 10.0)
    triangle = [
        SVector(1.0, 1.0, 1.0),
        SVector(2.0, 1.0, 1.0),
        SVector(1.0, 2.0, 1.0),
    ]
    equal_masses = ones(3)
    unequal_masses = [1.0, 2.0, 1.0]

    @test molecule_center_of_mass(triangle, equal_masses, box) ≈
          SVector(4 / 3, 4 / 3, 1.0)
    @test molecule_center_of_mass(triangle, unequal_masses, box) ≈
          SVector(1.5, 1.25, 1.0)
    @test orientation(CenterToAtomOrientation(1), triangle, equal_masses, box) ≈
          normalize(SVector(-1.0, -1.0, 0.0))
    @test orientation(CenterToAtomOrientation(1), triangle, unequal_masses, box) ≈
          normalize(SVector(-0.5, -0.25, 0.0))
    @test orientation(PlaneNormalOrientation(), triangle, box) ≈
          SVector(0.0, 0.0, 1.0)
    @test orientation(PlaneNormalOrientation(1, 3, 2), triangle, box) ≈
          SVector(0.0, 0.0, -1.0)

    wrapped_triangle = [
        SVector(9.5, 1.0, 1.0),
        SVector(0.5, 1.0, 1.0),
        SVector(9.5, 2.0, 1.0),
    ]
    @test orientation(CenterToAtomOrientation(1), wrapped_triangle, equal_masses, box) ≈
          normalize(SVector(-1.0, -1.0, 0.0))
    @test orientation(PlaneNormalOrientation(), wrapped_triangle, box) ≈
          SVector(0.0, 0.0, 1.0)

    @test_throws ArgumentError PlaneNormalOrientation(1, 1, 3)
    @test_throws ArgumentError molecule_center_of_mass(triangle, [1.0, 2.0], box)
    @test_throws ArgumentError molecule_center_of_mass(triangle, [1.0, 0.0, 1.0], box)
    @test_throws ArgumentError orientation(
        PlaneNormalOrientation(),
        fill(SVector(1.0, 1.0, 1.0), 3),
        box,
    )
end

@testset "Aligning field with displacement" begin
    common_args = Dict(
        "temperature" => [2.0],
        "model" => ["Trimer"],
        "list_type" => "LinkedList",
    )
    field = AligningField(SVector(0.0, 0.0, 0.7), PlaneNormalOrientation())
    no_field_system = load_chains("molecule.exyz", args=common_args)[1]
    field_system = load_chains(
        "molecule.exyz",
        args=merge(common_args, Dict("external_field" => field)),
    )[1]

    @test field_system.energy[1] - no_field_system.energy[1] ≈
          total_field_energy(field, field_system)

    molecule_id = field_system.molecule[1]
    old_field_energy = field_energy(field, field_system, molecule_id)
    no_field_action = Displacement(1, SVector(0.0, 0.0, 0.02), 0.0)
    field_action = Displacement(1, SVector(0.0, 0.0, 0.02), 0.0)
    ParticlesMC.perform_action!(no_field_system, no_field_action)
    ParticlesMC.perform_action!(field_system, field_action)
    new_field_energy = field_energy(field, field_system, molecule_id)

    @test field_action.δe - no_field_action.δe ≈
          new_field_energy - old_field_energy atol=1e-10
    @test field_system.external_field === field

    config_field = ParticlesMC.external_field_from_config(Dict(
        "type" => "align",
        "h" => [0.0, 0.0, 0.7],
        "orientation" => "center_to_atom",
        "atom" => 2,
    ))
    @test config_field isa AligningField
    @test config_field.orientation_definition == CenterToAtomOrientation(2)
end

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
    chains_vl = load_chains("config_0.lmp", args=Dict("temperature" => [0.231], "model" => ["JBB"], "list_type" => "VerletList", "list_parameters" => Dict("dr" => 0.2)))
    system_el = chains_el[1]
    system_ll = chains_ll[1]
    system_vl = chains_vl[1]
    @test system_el.N == system_ll.N
    @test system_el.d == system_ll.d
    @test system_el.temperature == system_ll.temperature
    @test all.(system_el.position == system_ll.position)
    @test all.(system_el.species == system_ll.species)
    energy_el = system_el.energy[1] / length(system_el)
    energy_ll = system_ll.energy[1] / length(system_ll)
    energy_vl = system_vl.energy[1] / length(system_vl)
    @test isapprox(energy_el, -2.676832, atol=1e-6)
    @test isapprox(energy_ll, -2.676832, atol=1e-6)
    @test isapprox(energy_vl, -2.676832, atol=1e-6)

    # Test simulation energy
    M = 1
    seed = 10
    sp1, sp2, sp3 = findall(isequal(1), system_el.species), findall(isequal(2), system_el.species), findall(isequal(3), system_el.species)
    NA, NB, NC = length(sp1), length(sp2), length(sp3)
    steps = 100
    burn = 0
    block = [0, 1, 2, 4, 8]
    sampletimes = build_schedule(steps, burn, block)

    # NO SWAPS
    pswap = 0.0
    displacement_policy = SimpleGaussian()
    displacement_parameters = ComponentArray(σ=0.05)
    pool = (
        Move(Displacement(0, zero(system_el.box), 0.0), displacement_policy, displacement_parameters, 1 - pswap),
    )
    algorithm_list = (
    (algorithm=Metropolis, pool=pool, seed=seed, parallel=false, sweepstep=system_el.N),
    (algorithm=StoreCallbacks, callbacks=(energy,), scheduler=sampletimes),
    (algorithm=StoreAcceptance, dependencies=(Metropolis,), scheduler=sampletimes), 
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

    ## Verlet List simulation
    chains_vl = [deepcopy(system_vl)]
    path_vl = "data/noswap/verlet_list/"
    simulation = Simulation(chains_vl, algorithm_list, steps; path=path_vl, verbose=true)
    run!(simulation)

    ## Read energy data and compare
    path_energy_el = joinpath(path_el, "chains/1/energy.dat")
    path_energy_ll = joinpath(path_ll, "chains/1/energy.dat")
    path_energy_vl = joinpath(path_vl, "chains/1/energy.dat")
    energy_el= readdlm(path_energy_el)[:, 2]
    energy_ll = readdlm(path_energy_ll)[:, 2]
    energy_vl = readdlm(path_energy_vl)[:, 2]
    @test isapprox(energy_el, energy_ll, atol=1e-6)
    @test isapprox(energy_ll, energy_vl, atol=1e-6)

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
        (algorithm=StoreCallbacks, callbacks=(energy,), scheduler=sampletimes),
        (algorithm=StoreAcceptance, dependencies=(Metropolis,), scheduler=sampletimes),
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
    path_energy_el = joinpath(path_el, "chains/1/energy.dat")
    path_energy_ll = joinpath(path_ll, "chains/1/energy.dat")
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

    # NO SWAPS
    pswap = 0.0
    displacement_policy = SimpleGaussian()
    displacement_parameters = ComponentArray(σ=0.05)
    pool = (
        Move(Displacement(0, zero(system_el.box), 0.0), displacement_policy, displacement_parameters, 1 - pswap),
    )
    algorithm_list = (
    (algorithm=Metropolis, pool=pool, seed=seed, parallel=false, sweepstep=system_el.N),
    (algorithm=StoreCallbacks, callbacks=(energy,), scheduler=sampletimes),
    (algorithm=StoreAcceptance, dependencies=(Metropolis,), scheduler=sampletimes),        
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
    path_energy_el = joinpath(path_el, "chains/1/energy.dat")
    path_energy_ll = joinpath(path_ll, "chains/1/energy.dat")
    energy_el= readdlm(path_energy_el)[:, 2]
    energy_ll = readdlm(path_energy_ll)[:, 2]
    @test isapprox(energy_el, energy_ll, atol=1e-6)

end

@testset "molecules rotation test" begin
    chains = load_chains("molecule.exyz", args=Dict("temperature" => [2.0], "model" => ["Trimer"], "list_type" => "LinkedList"))
    system = chains[1]
    steps = 128
    burn  = 0
    theta_T = [0.0, π/4, π]
    sampletimes = build_schedule(steps, burn, 2.0)
    displacement_policy     = SimpleGaussian()
    displacement_parameters = ComponentArray(σ=0.05)
    pool = (Move(Displacement(0, zero(system.box), 0.0), displacement_policy, displacement_parameters, 1.0),)
    algorithm_list = (
        (algorithm=Metropolis, pool=pool, seed=10, parallel=false, sweepstep=system.N),
        (algorithm=ComputeRotation, scheduler=build_schedule(steps, burn, 1), theta_T=theta_T),
        (algorithm=StorePhiTrajectories, scheduler=sampletimes, path="data/rotation/"),
        (algorithm=StoreLastPhiFrame, scheduler=[steps], path="data/rotation/"),
    )
    chains = [deepcopy(system)]
    simulation = Simulation(chains, algorithm_list, steps; path="data/rotation/", verbose=true)
    run!(simulation)
    println("OK")
end

@testset "molecules rotation test" begin
    chains = load_chains("molecule.exyz", args=Dict("temperature" => [2.0], "model" => ["Trimer"], "list_type" => "LinkedList"))
    system = chains[1]
    steps = 128
    burn  = 0
    theta_T = [0.0, π/4, π]
    sampletimes = build_schedule(steps, burn, 2.0)
    displacement_policy     = SimpleGaussian()
    displacement_parameters = ComponentArray(σ=0.05)
    pool = (Move(Displacement(0, zero(system.box), 0.0), displacement_policy, displacement_parameters, 1.0),)
    algorithm_list = (
        (algorithm=Metropolis, pool=pool, seed=10, parallel=false, sweepstep=system.N),
        (algorithm=ComputeRotation, scheduler=build_schedule(steps, burn, 1), theta_T=theta_T),
        (algorithm=StorePhiTrajectories, scheduler=sampletimes, path="data/rotation/"),
        (algorithm=StoreLastPhiFrame, scheduler=[steps], path="data/rotation/"),
    )
    chains = [deepcopy(system)]
    simulation = Simulation(chains, algorithm_list, steps; path="data/rotation/", verbose=true)
    run!(simulation)
    println("OK")
end
