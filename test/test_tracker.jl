using ParticlesMC, StaticArrays, LinearAlgebra

# check MSADState constructs empty
s = ParticlesMC.MSADState{Float64}()
println("MSADState empty: ", s.initialized == false)

# check MSADTracker constructor with valid schedules
compute = collect(0:10:100)
output  = collect(0:20:100)
chains  = [1, 2]
t = ParticlesMC.MSADTracker(chains, pi/4, compute, output, "/tmp/test_msad")
println("MSADTracker created: ", t.theta_T ≈ pi/4)

# check assertion fires on invalid schedules
try
    bad_output = [0, 5, 10]   # 5 not in compute schedule
    ParticlesMC.MSADTracker(chains, pi/4, compute, bad_output, "/tmp/test_msad")
    println("ERROR: should have thrown")
catch e
    println("Assertion caught correctly: ", e isa AssertionError)
end