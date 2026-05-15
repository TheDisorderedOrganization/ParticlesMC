using StaticArrays, LinearAlgebra, BenchmarkTools

function rotation_vector_skew(R::SMatrix{3,3,T}) where {T}
    cos_θ = clamp((tr(R) - 1) / 2, -one(T), one(T))
    skew  = (R - R') / 2
    sin_n = SVector(skew[3,2], skew[1,3], skew[2,1])
    sin_θ = norm(sin_n)
    θ     = atan(sin_θ, cos_θ)
    if sin_θ > 1e-6
        return (θ / sin_θ) * sin_n
    elseif cos_θ > 0
        return sin_n
    else
        RpI = R + one(R)
        c1  = SVector(RpI[1,1], RpI[2,1], RpI[3,1])
        c2  = SVector(RpI[1,2], RpI[2,2], RpI[3,2])
        c3  = SVector(RpI[1,3], RpI[2,3], RpI[3,3])
        n   = norm(c1) >= norm(c2) ? c1 : c2
        n   = norm(n)  >= norm(c3) ? n  : c3
        return π * (n / norm(n))
    end
end

function rotation_vector_eig(R::SMatrix{3,3,T}) where {T}
    cos_θ = clamp((tr(R) - 1) / 2, -one(T), one(T))
    vals, vecs = eigen(R)
    idx = argmin(abs.(real.(vals) .- 1))
    n = real.(vecs[:, idx])
    n = SVector{3,T}(n[1], n[2], n[3])
    n_skew = SMatrix{3,3,T}(
         0,     n[3], -n[2],
        -n[3],  0,     n[1],
         n[2], -n[1],  0
    )
    sin_θ = clamp(-tr(n_skew * R) / 2, -one(T), one(T))
    θ = atan(sin_θ, cos_θ)
    return θ * n
end

# 90° rotation around z — normal case
R = SMatrix{3,3,Float64}(0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0)

println("=== Skew-symmetric ===")
@btime rotation_vector_skew($R)

println("=== Eigenvector ===")
@btime rotation_vector_eig($R)
