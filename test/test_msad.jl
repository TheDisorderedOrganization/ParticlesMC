using ParticlesMC
using StaticArrays
using LinearAlgebra
using Test

@testset "body_frame" begin

    # known inputs — worked out by hand
    r1 = SVector(0.0, 0.0, 0.0)
    r2 = SVector(1.0, 0.0, 0.0)
    r3 = SVector(0.0, 1.0, 0.0)
    L  = 10.0

    R = ParticlesMC.body_frame(r1, r2, r3, L)

    # expected result from hand calculation:
    # e1 = [1,  0,  0]
    # e2 = [0,  0,  1]
    # e3 = [0, -1,  0]
    @test R[:, 1] ≈ SVector( 1.0,  0.0,  0.0)
    @test R[:, 2] ≈ SVector( 0.0,  0.0,  1.0)
    @test R[:, 3] ≈ SVector( 0.0, -1.0,  0.0)

    # R must be a valid rotation matrix:
    # columns orthonormal → R'R = I
    # right-handed → det(R) = +1
    @test R' * R ≈ I
    @test det(R) ≈ 1.0

end

@testset "body_frame minimum image" begin

    # same molecule but r2 is wrapped across the boundary
    # r1 at 0.1, r2 at 9.9 — bond should be -0.2, not 9.8
    r1 = SVector(0.1, 0.0, 0.0)
    r2 = SVector(9.9, 0.0, 0.0)
    r3 = SVector(0.1, 1.0, 0.0)
    L  = 10.0

    R_wrapped = ParticlesMC.body_frame(r1, r2, r3, L)

    # unwrapped version — same molecule, same frame
    r2_unwrapped = SVector(-0.1, 0.0, 0.0)
    R_unwrapped  = ParticlesMC.body_frame(r1, r2_unwrapped, r3, L)

    @test R_wrapped ≈ R_unwrapped

end

@testset "rotation_vector" begin

    # identity matrix → zero rotation
    R_id = SMatrix{3,3}(1.0I)
    @test ParticlesMC.rotation_vector(R_id) ≈ SVector(0.0, 0.0, 0.0)

    # 90° rotation around z-axis
    # R = [0 -1 0]
    #     [1  0 0]
    #     [0  0 1]
    R_90z = SMatrix{3,3}(0.0, 1.0, 0.0,
                        -1.0, 0.0, 0.0,
                         0.0, 0.0, 1.0)
    Ω = ParticlesMC.rotation_vector(R_90z)
    # should be θ=π/2 around z → Ω = [0, 0, π/2]
    @test norm(Ω) ≈ π/2
    @test Ω / norm(Ω) ≈ SVector(0.0, 0.0, 1.0)

end

@testset "body_frame PBC - molecule split across each axis" begin

    L = 10.0

    # --- x axis ---
    # r1 near left wall, r2 wrapped to right wall
    # true bond is [-0.2, 0, 0], not [9.8, 0, 0]
    r1_x = SVector(0.1, 0.0, 0.0)
    r2_x = SVector(9.9, 0.0, 0.0)   # wrapped
    r3_x = SVector(0.1, 1.0, 0.0)
    R_wrapped_x   = ParticlesMC.body_frame(r1_x, r2_x, r3_x, L)
    R_unwrapped_x = ParticlesMC.body_frame(r1_x, SVector(-0.1, 0.0, 0.0), r3_x, L)
    @test R_wrapped_x ≈ R_unwrapped_x

    # --- y axis ---
    r1_y = SVector(0.0, 0.1, 0.0)
    r2_y = SVector(1.0, 0.1, 0.0)
    r3_y = SVector(0.0, 9.9, 0.0)   # r3 wrapped across y
    R_wrapped_y   = ParticlesMC.body_frame(r1_y, r2_y, r3_y, L)
    R_unwrapped_y = ParticlesMC.body_frame(r1_y, r2_y, SVector(0.0, -0.1, 0.0), L)
    @test R_wrapped_y ≈ R_unwrapped_y

    # --- z axis ---
    r1_z = SVector(0.0, 0.0, 0.1)
    r2_z = SVector(1.0, 0.0, 0.1)
    r3_z = SVector(0.0, 1.0, 0.1)
    r2_z_wrapped = SVector(1.0, 0.0, 9.9 + 0.1)   # wrapped across z... wait
    # r1 near top wall, r2 wrapped to bottom
    r1_z2 = SVector(0.0, 0.0, 0.1)
    r2_z2 = SVector(0.0, 0.0, 9.9)   # wrapped
    r3_z2 = SVector(1.0, 0.0, 0.1)
    R_wrapped_z   = ParticlesMC.body_frame(r1_z2, r2_z2, r3_z2, L)
    R_unwrapped_z = ParticlesMC.body_frame(r1_z2, SVector(0.0, 0.0, -0.1), r3_z2, L)
    @test R_wrapped_z ≈ R_unwrapped_z

end

@testset "body_frame PBC - result is always a valid rotation matrix" begin

    L = 10.0
    # test many random molecules, some crossing boundaries
    # a valid rotation matrix always satisfies R'R = I and det(R) = +1
    import Random
    Random.seed!(42)
    for _ in 1:100
        # random positions anywhere in box
        r1 = SVector(rand()*L, rand()*L, rand()*L)
        r2 = SVector(rand()*L, rand()*L, rand()*L)
        r3 = SVector(rand()*L, rand()*L, rand()*L)
        R = ParticlesMC.body_frame(r1, r2, r3, L)
        @test R' * R ≈ I  atol=1e-10
        @test det(R) ≈ 1.0  atol=1e-10
    end

end

@testset "rotation_vector PBC - roundtrip" begin

    # if we compute R between two frames, rotation_vector should
    # give back an Ω whose norm equals the angle between them
    # and applying it back should reconstruct the rotation

    # 180° rotation around x axis
    R_180x = SMatrix{3,3}(1.0,  0.0,  0.0,
                          0.0, -1.0,  0.0,
                          0.0,  0.0, -1.0)
    Ω = ParticlesMC.rotation_vector(R_180x)
    @test norm(Ω) ≈ π  atol=1e-10
    @test Ω / norm(Ω) ≈ SVector(1.0, 0.0, 0.0)  atol=1e-10

    # 180° rotation around y axis
    R_180y = SMatrix{3,3}(-1.0, 0.0, 0.0,
                           0.0, 1.0, 0.0,
                           0.0, 0.0,-1.0)
    Ω = ParticlesMC.rotation_vector(R_180y)
    @test norm(Ω) ≈ π  atol=1e-10
    @test Ω / norm(Ω) ≈ SVector(0.0, 1.0, 0.0)  atol=1e-10

    # small rotation — near-identity, tests the sin≈0 branch
    ε = 1e-8
    R_tiny = SMatrix{3,3}(1.0,  ε,   0.0,
                          -ε,  1.0,  0.0,
                          0.0, 0.0,  1.0)
    Ω = ParticlesMC.rotation_vector(R_tiny)
    @test norm(Ω) < 1e-6   # nearly zero rotation

end