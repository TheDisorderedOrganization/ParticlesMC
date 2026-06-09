using LinearAlgebra, StaticArrays, BenchmarkTools

# ── original ──────────────────────────────────────────────────────────────────
function rotation_vector_eigen(R::SMatrix{3,3,T}) where {T}
    cos_θ = clamp((tr(R) - 1) / 2, -one(T), one(T))
    vals, vecs = eigen(R)
    idx = argmin(abs.(real.(vals) .- 1))
    n = real.(vecs[:, idx])
    n = SVector{3,T}(n[1], n[2], n[3])
    n_skew = SMatrix{3,3,T}(0, n[3], -n[2], -n[3], 0, n[1], n[2], -n[1], 0)
    sin_θ = clamp(-tr(n_skew * R) / 2, -one(T), one(T))
    θ = atan(sin_θ, cos_θ)
    return θ * n
end

# ── closed-form ───────────────────────────────────────────────────────────────
function rotation_vector_closed(R::SMatrix{3,3,T}) where {T}
    cos_θ = clamp((tr(R) - 1) / 2, -one(T), one(T))
    θ = acos(cos_θ)

    ax = R[3,2] - R[2,3]
    ay = R[1,3] - R[3,1]
    az = R[2,1] - R[1,2]
    sin_θ = sqrt(ax^2 + ay^2 + az^2) / 2   # exact: norm of 2·sin(θ)·n / 2

    if θ < sqrt(eps(T))                      # θ ≈ 0 : no rotation
        return zero(SVector{3,T})
    elseif π - θ < sqrt(eps(T))             # θ ≈ π : antisymmetric part vanishes
        # use diagonal of symmetric part: (R+I)/2 = n⊗n
        nx = sqrt(max((R[1,1] + 1) / 2, zero(T)))
        ny = sqrt(max((R[2,2] + 1) / 2, zero(T)))
        nz = sqrt(max((R[3,3] + 1) / 2, zero(T)))
        n  = SVector{3,T}(nx, ny, nz)
        return θ * n / norm(n)
    else                                     # general case: zero allocs
        inv_2s = θ / (2 * sin_θ)
        return SVector{3,T}(ax * inv_2s, ay * inv_2s, az * inv_2s)
    end
end

# ── helper: build exact rotation matrix from axis + angle ─────────────────────
function make_rotation(axis::SVector{3,Float64}, θ::Float64)
    n = axis / norm(axis)
    # skew matrix stored column-major for SMatrix
    N = SMatrix{3,3,Float64}(
         0,    n[3], -n[2],
        -n[3],  0,    n[1],
         n[2], -n[1],  0
    )
    return SMatrix{3,3,Float64}(I + sin(θ)*N + (1 - cos(θ))*(N*N))
end

# ── Test 1: 10_000 random rotations, θ ∈ (0.01, π-0.01) ─────────────────────
println("="^60)
println("Test 1 — 10_000 random rotations, θ ∈ (0.01, π-0.01)")
println("="^60)
max_err = 0.0
for _ in 1:10_000
    axis   = SVector{3,Float64}(randn(), randn(), randn())
    θ_true = 0.01 + (π - 0.02) * rand()
    R      = make_rotation(axis, θ_true)
    v1     = rotation_vector_eigen(R)
    v2     = rotation_vector_closed(R)
    err    = norm(v1 - v2)
    global max_err = max(max_err, err)
end
println("Max |eigen - closed|: ", max_err, "  (should be < 1e-12)")

# ── Test 2: near-zero (θ → 0) ────────────────────────────────────────────────
println()
println("="^60)
println("Test 2 — near-zero rotation θ ≈ 0")
println("="^60)
for θ_small in [1e-4, 1e-6, 1e-8, 1e-12]
    axis = SVector{3,Float64}(1.0, 0.0, 0.0)
    R    = make_rotation(axis, θ_small)
    v1   = rotation_vector_eigen(R)
    v2   = rotation_vector_closed(R)
    println("θ=%.0e  |v_eigen|=%.6e  |v_closed|=%.6e  diff=%.2e\n",
            θ_small, norm(v1), norm(v2), norm(v1-v2))
end

# ── Test 3: near-π — RELEVANT because theta_T=3.0 rad ────────────────────────
println()
println("="^60)
println("Test 3 — near-π rotation (theta_T=3.0 means you DO reach this)")
println("="^60)
for δ in [1e-1, 1e-3, 1e-5, 1e-7]
    axis = normalize(SVector{3,Float64}(1.0, 2.0, 3.0))
    R    = make_rotation(axis, π - δ)
    v1   = rotation_vector_eigen(R)
    v2   = rotation_vector_closed(R)
    println("θ=π-%.0e  |v_eigen|=%.6f  |v_closed|=%.6f  diff=%.2e\n",
            δ, norm(v1), norm(v2), norm(v1-v2))
end

# ── Test 4: benchmark ─────────────────────────────────────────────────────────
println()
println("="^60)
println("Test 4 — benchmark: 1000 calls each")
println("="^60)
Rs = [make_rotation(normalize(SVector{3,Float64}(randn(),randn(),randn())),
                    0.1 + (π-0.2)*rand())
      for _ in 1:1000]

t_eigen  = @belapsed sum(norm, rotation_vector_eigen.($(Rs)))
t_closed = @belapsed sum(norm, rotation_vector_closed.($(Rs)))
println("eigen  : %8.2f µs\n", t_eigen  * 1e6)
println("closed : %8.2f µs\n", t_closed * 1e6)
println("speedup: %.1fx\n", t_eigen / t_closed)

# ── Test 5: verify the identity R - Rᵀ = 2 sin(θ) N analytically ─────────────
println()
println("="^60)
println("Test 5 — verify R - Rᵀ = 2·sin(θ)·N identity directly")
println("="^60)
axis  = normalize(SVector{3,Float64}(1.0, 2.0, 3.0))
θ_ref = 1.23456
R     = make_rotation(axis, θ_ref)
antisym = (R - R') / 2                       # = sin(θ)·N
n_from_antisym = SVector{3,Float64}(
    antisym[3,2],   # = sin(θ)·n[1]
    antisym[1,3],   # = sin(θ)·n[2]
    antisym[2,1],   # = sin(θ)·n[3]
)
n_recovered = n_from_antisym / norm(n_from_antisym)
println("True axis     : ", round.(axis,          sigdigits=8))
println("Recovered axis: ", round.(n_recovered,   sigdigits=8))
println("Error on axis : ", norm(axis - n_recovered))
println("True θ        : ", θ_ref)
println("Recovered θ   : ", atan(norm(n_from_antisym), (tr(R)-1)/2))