using MonteCarlo
using ParticlesMC
using Distributions
using Random
using StaticArrays
using ComponentArrays

seed = 42
rng = Xoshiro(seed)
N = 3000
M = 1
Length = 3
d = 3
temperature = 2.0
density = 1.2
box = @SVector fill(typeof(temperature)((N / density)^(1 / d)), d)
position_1 = [box .* @SVector rand(rng, d) for i in 1:Int(N // Length), m in 1:M]
position = []
for r in position_1
    push!(position, r)
    push!(position, r .+ @SVector [0.1, 0.1, 0.1])
    push!(position, r .- @SVector [0.1, 0.1, 0.1])
end
position = Vector{SVector{3,Float64}}(position)
mol_species = Vector{Tuple{Int, Int, Int, Int}}()
species = Vector{Int}()
for i in 1:Int(N // Length)
    push!(species, 1)
    push!(species, 2)
    push!(species, 3)
    push!(mol_species, (1, 3, 3 * (i - 1) + 1, 3 * i))
end

function create_bond_matrix(N::Int)
    # Create a vector to store the SVectors, each containing a pair of integers
    matrix = Vector{Vector{Int}}()
    # Populate the matrix with pairs according to the specified pattern
    for i in 2:3:N
        push!(matrix, [i, i + 1])
        push!(matrix, [i - 1, i + 1])
        push!(matrix, [i - 1, i])
    end
    return matrix
end

bonds = create_bond_matrix(N)
epsilon = SMatrix{3, 3, Float64}([1.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 1.0])
sigma = SMatrix{3, 3, Float64}([0.9 0.95 1.0; 0.95 1.0 1.05; 1.0 1.05 1.1])
k = SMatrix{3, 3, Float64}([0.0 33.241 30.0; 33.241 0.0 27.210884; 30.0 27.210884 0.0])
r0 = SMatrix{3, 3, Float64}([0.0 1.425 1.5; 1.425 0.0 1.575; 1.5 1.575 0.0])
model = GeneralKG(epsilon, sigma, k, r0)
chains = [System(position, species, mol_species, density, temperature, model, bonds; list_type=LinkedList) for _ in 1:M]
## Define moves and combine them into M independent pools
pswap = 0.2
displacement_policy = SimpleGaussian()
displacement_parameters = ComponentArray(σ=0.05)
pools = [(
    Move(Displacement(0, zero(box)), displacement_policy, displacement_parameters, 1.0),
) for _ in 1:M]
## Define the simulation struct
steps = 1000
burn = 0
block = [0, 10]
sampletimes = build_schedule(steps, burn, block)
schedulers = [build_schedule(steps, 0, 1), sampletimes, sampletimes, [0, steps], build_schedule(steps, burn, steps ÷ 10)]
callbacks = (callback_energy, callback_acceptance)

path = "data/test/particles/Molecules/T$temperature/N$N/M$M/seed$seed"
algorithms = (
    Metropolis(chains, pools; sweepstep=N, seed=seed, parallel=false),
    StoreCallbacks(callbacks, path),
    StoreTrajectories(chains, path),
    StoreLastFrames(chains, path),
    PrintTimeSteps(),
    )
## Run the simulation :)
simulation = Simulation(chains, algorithms, steps; schedulers=schedulers, path=path, verbose=true)
run!(simulation)