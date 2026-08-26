using Arianna
using ComponentArrays
using LinearAlgebra
using ParticlesMC
using Statistics
using StaticArrays

const FIELD_DIRECTION = SVector(0.0, 0.0, 1.0)

# This callback deliberately measures molecule 1: the example contains exactly
# one molecule, so no molecule-selection feature is needed for this first test.
Arianna.@callback function field_alignment(system::Molecules)
    n = orientation(PlaneNormalOrientation(), system, 1)
    return dot(FIELD_DIRECTION, n)
end

function isolated_otp(field_magnitude)
    # The initial triangle lies in the yz plane, so its normal starts
    # perpendicular to the field direction (+z).
    positions = [
        SVector(2.0, 2.0, 2.0),
        SVector(2.0, 2.8, 2.0),
        SVector(2.0, 2.4, 2.7),
    ]
    species = [1, 2, 3]
    molecule = ones(Int, 3)
    bonds = [[2, 3], [1, 3], [1, 2]]
    box_length = 5.0
    density = length(positions) / box_length^3
    temperature = 1.0
    field = AligningField(
        field_magnitude * FIELD_DIRECTION,
        PlaneNormalOrientation(),
    )

    return System(
        positions,
        species,
        molecule,
        density,
        temperature,
        Trimer(),
        bonds;
        masses=ones(3),
        external_field=field,
        list_type=EmptyList,
    )
end

function read_callback_values(path)
    lines = filter(!isempty, readlines(path))
    return [parse(Float64, split(line)[2]) for line in lines]
end

function run_one_field(field_magnitude;
                       steps=50_000,
                       burn=10_000,
                       interval=20,
                       trajectory_interval=100)
    system = isolated_otp(field_magnitude)
    policy = SimpleGaussian()
    parameters = ComponentArray(σ=0.05)
    pool = (
        Move(Displacement(0, zero(system.box), 0.0), policy, parameters, 1.0),
    )
    scheduler = build_schedule(steps, burn, interval)
    trajectory_scheduler = build_schedule(steps, 0, trajectory_interval)
    output = joinpath(@__DIR__, "output", "h_$(field_magnitude)")
    algorithms = (
        (algorithm=Metropolis, pool=pool, seed=1234, parallel=false, sweepstep=system.N),
        (algorithm=StoreCallbacks, callbacks=(field_alignment,), scheduler=scheduler),
        (algorithm=StoreTrajectories, scheduler=trajectory_scheduler, fmt=XYZ()),
    )
    simulation = Simulation([system], algorithms, steps; path=output, verbose=false)
    run!(simulation)

    values = read_callback_values(joinpath(output, "chains", "1", "field_alignment.dat"))
    return mean(values), std(values) / sqrt(length(values)), length(values)
end

field_magnitudes = [0.0, 0.5, 1.0, 3.0, 10.0]
results = [run_one_field(h) for h in field_magnitudes]

summary_path = joinpath(@__DIR__, "alignment-vs-field.csv")
open(summary_path, "w") do io
    println(io, "h,mean_alignment,standard_error,n_samples")
    for (h, (mean_alignment, standard_error, n_samples)) in zip(field_magnitudes, results)
        println(io, "$h,$mean_alignment,$standard_error,$n_samples")
        println("h=$h: <hhat dot n> = $mean_alignment +/- $standard_error")
    end
end
println("Wrote $summary_path")
