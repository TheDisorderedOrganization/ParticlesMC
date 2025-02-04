using MonteCarlo
using ParticlesMC
using Test
using StaticArrays
using Distributions
using ComponentArrays
using DelimitedFiles

@testset "Potential energy test" begin
    # Test inital configuration
    data = readlines("test/config_0.xyz")
    N = parse(Int, data[4])
    gbox = parse.(Float64, split(data[6], " "))
    L = gbox[2] - gbox[1]
    box = @SVector [L, L]
    frame = data[10:end]
    species = zeros(Int, N)
    position = Vector{SVector{2,Float64}}(undef, N)
    for i in eachindex(frame)
        species[i] = parse(Int, split(frame[i], " ")[1])
        position[i] = fold_back(SVector{2,Float64}(parse.(Float64, split(frame[i], " ")[2:3])), box)
    end
    temperature = 0.231
    density = N / L^2
    system_el = System(position, species, density, temperature, JBB())
    system_ll = System(position, species, density, temperature, JBB(); list_type=LinkedList)
    energy_el = mean(system_el.local_energy) / 2
    energy_ll = mean(system_ll.local_energy) / 2
    @test isapprox(energy_el, -2.676832, atol=1e-6)
    @test isapprox(energy_ll, -2.676832, atol=1e-6)

    # Test simulation energy
    M = 1
    seed = 10
    sp1, sp2, sp3 = findall(isequal(1), species), findall(isequal(2), species), findall(isequal(3), species)
    NA, NB, NC = length(sp1), length(sp2), length(sp3)
    steps = 100
    burn = 0
    block = [0, 1, 2, 4, 8]
    sampletimes = scheduler(steps, burn, block)
    callbacks = (callback_energy, callback_acceptance)

    # NO SWAPS
    pswap = 0.0
    displacement_policy = SimpleGaussian()
    displacement_parameters = ComponentArray(σ=0.05)
    pools = [(
        Move(Displacement(0, zero(box)), displacement_policy, displacement_parameters, 1 - pswap),
    ) for _ in 1:M]
    
    ## Empty List simulation
    chains_el = [deepcopy(system_el)]
    path_el = "data/noswap/empty_list/"
    simulation = Simulation(chains_el, pools, steps; sweepstep=N, sampletimes=sampletimes, seed=seed, parallel=false, verbose=true, path=path_el)
    run!(simulation, callbacks...)
    
    ## Linked List simulation
    chains_ll = [deepcopy(system_ll)]
    path_ll = "data/noswap/linked_list/"
    simulation = Simulation(chains_ll, pools, steps; sweepstep=N, sampletimes=sampletimes, seed=seed, parallel=false, verbose=true, path=path_ll)
    run!(simulation, callbacks...)

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
    pools = [(
    Move(Displacement(0, zero(box)), displacement_policy, displacement_parameters, 1 - pswap),
    Move(DiscreteSwap(0, 0, (1, 3), (NA, NC)), swap_policy, swap_parameters, pswap / 2),
    Move(DiscreteSwap(0, 0, (2, 3), (NB, NC)), swap_policy, swap_parameters, pswap / 2),
    ) for _ in 1:M]

    ## Empty List simulation
    chains_el = [deepcopy(system_el)]
    path_el = "data/swap/empty_list/"
    simulation = Simulation(chains_el, pools, steps; sweepstep=N, sampletimes=sampletimes, seed=seed, parallel=false, verbose=true, path=path_el)
    run!(simulation, callbacks...)
    
    ## Linked List simulation
    chains_ll = [deepcopy(system_ll)]
    path_ll = "data/swap/linked_list/"
    simulation = Simulation(chains_ll, pools, steps; sweepstep=N, sampletimes=sampletimes, seed=seed, parallel=false, verbose=true, path=path_ll)
    run!(simulation, callbacks...)

    ## Read energy data and compare
    path_energy_el = joinpath(path_el, "energy.dat")
    path_energy_ll = joinpath(path_ll, "energy.dat")
    energy_el= readdlm(path_energy_el)[:, 2]
    energy_ll = readdlm(path_energy_ll)[:, 2]
    @test isapprox(energy_el, energy_ll, atol=1e-6)


end