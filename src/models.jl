using LinearAlgebra, StaticArrays

"""
    abstract type Model end

An abstract base type representing a general interaction model.
"""
abstract type Model end

"""
    abstract type DiscreteModel end

An abstract base type representing a general interaction model with a fixed number of species.
"""
abstract type DiscreteModel <: Model end


"""
    abstract type ContinuousModel end

An abstract base type representing a general interaction model with a continuous size polidispersity.
"""
abstract type ContinuousModel <: Model end

## POTENTIALS
harmonic_spheres(r2, ϵ, σ2) = ϵ * (1 - r2 / σ2)

inverse_power(r2, ϵ, σ2, ndiv2) = ϵ * (σ2 / r2)^(ndiv2)

@inline function lennard_jones(r2::T, ϵ4::T, σ2::T) where T<:AbstractFloat
    x = σ2 * inv(r2)
    x3 = x * x * x  # (σ²/r²)^3
    return ϵ4 * (x3 * x3 - x3)
end

fene(r2, kr02, r02) = kr02 * log(1 - r2 * inv(r02))

###############################################################################
"""
    struct SoftSpheres{T<:AbstractArray} <: DiscreteModel

A struct representing a soft-sphere interaction model where particles interact via an inverse power-law potential.

# Fields
- `name::String`: A descriptive name for the model.
- `ϵ::T`: Energy scale matrix for particle interactions between different types.
- `σ::T`: Characteristic length matrix for particle interactions between different types.
- `n::Int`: Exponent controlling the softness of the interaction.
- `rcut::T`: Cutoff distance matrix beyond which interactions are neglected.
- `shift::T`: Shift matrix to ensure the potential is zero at the cutoff distance.
"""
struct SoftSpheres{T<:AbstractArray, N<:Number} <: DiscreteModel
    name::String
    ϵ::T
    σ::T
    σ2::T
    n::Int
    ndiv2::N
    rcut::T
    rcut2::T
    shift::T
end

function SoftSpheres(ϵ, σ, n; rcut=2.5*σ, name="SoftSpheres")
    σ2 = σ ^ 2
    rcut2 = rcut ^ 2
    ndiv2 = isodd(n) ? n / 2 : n ÷ 2
    shift = inverse_power(rcut2, ϵ, σ2, ndiv2)
    return SoftSpheres(name, ϵ, σ, σ2, n, ndiv2, rcut, rcut2, shift)
end

function potential(r2, model::SoftSpheres)
    return inverse_power(r2, model.ϵ, model.σ2, model.ndiv2) - model.shift
end

function BHHP()
    ϵ = [1.0 1.0; 1.0 1.0]
    σ = [1.0 1.2; 1.2 1.4]
    LJ_11 = SoftSpheres(ϵ[1,1], σ[1,1], 12)
    LJ_12 = SoftSpheres(ϵ[1,2], σ[1,2], 12)
    LJ_21 = SoftSpheres(ϵ[2,1], σ[2,1], 12)
    LJ_22 = SoftSpheres(ϵ[2,2], σ[2,2], 12)
    return [LJ_11 LJ_12; LJ_21 LJ_22]
end

###############################################################################
"""
    struct LennardJones{T<:AbstractArray} <: DiscreteModel

A struct representing the KobAndersen interaction model where particles interact via a Lennard-Jones potential.

# Fields
- `name::String`: A descriptive name for the model.
- `ϵ::T`: Energy scale matrix for particle interactions between different types.
- `σ::T`: Characteristic length matrix for particle interactions between different types.
- `rcut::T`: Cutoff distance matrix beyond which interactions are neglected.
- `shift::T`: Shift matrix to ensure the potential is zero at the cutoff distance.
"""
struct LennardJones{T<:AbstractFloat} <: DiscreteModel
    name::String
    ϵ::T
    ϵ4::T
    σ::T
    σ2::T
    rcut::T
    rcut2::T
    shift::T
end

function LennardJones(ϵ, σ; rcut=2.5*σ, name="LennardJones")
    σ2 = σ ^ 2
    rcut2 = rcut ^ 2
    shift = lennard_jones(rcut2, 4ϵ, σ2)
    return LennardJones(name, ϵ, 4ϵ, σ, σ2, rcut, rcut2, shift)
end

function potential(r2, model::LennardJones)
    return lennard_jones(r2, model.ϵ4, model.σ2) - model.shift
end

function KobAndersen()
    ϵ = [1.0 1.5; 1.5 0.5]
    σ = [1.0 0.8; 0.8 0.88]
    LJ_11 = LennardJones(ϵ[1,1], σ[1,1])
    LJ_12 = LennardJones(ϵ[1,2], σ[1,2])
    LJ_21 = LennardJones(ϵ[2,1], σ[2,1])
    LJ_22 = LennardJones(ϵ[2,2], σ[2,2])
    return [LJ_11 LJ_12; LJ_21 LJ_22]
end

###############################################################################

struct SmoothLennardJones{T<:AbstractFloat} <: DiscreteModel
    name::String
    ϵ::T
    σ::T
    ϵ4::T
    σ2::T
    C0::T
    C2_σ2::T
    C4_σ4::T
    rcut::T
    rcut2::T
end

function SmoothLennardJones(ϵ::T, σ::T; rcut::T=2.5*σ, name = "SmoothLennardJones") where T <: AbstractFloat
    C0 = T(0.04049023795)
    C2 = T(-0.00970155098)
    C4 = T(0.00062012616)
    σ2= σ ^ 2
    C2_σ2 = C2 / σ2
    C4_σ4 = C4 / σ2 ^ 2
    rcut2= rcut^2
    return SmoothLennardJones(name, ϵ, σ, 4ϵ, σ2, C0, C2_σ2, C4_σ4, rcut, rcut2)
end

function potential(r2::T, model::SmoothLennardJones) where T <: AbstractFloat
    ϵ4 = model.ϵ4
    lj = lennard_jones(r2, ϵ4, model.σ2)
    shift = ϵ4 * (model.C0 + r2 * muladd(r2, model.C4_σ4, model.C2_σ2))
    return lj + shift
end

function JBB()
    ϵ = [1.0 1.5 0.75; 1.5 0.5 1.5; 0.75 1.5 0.75]
    σ = [1.0 0.8 0.9; 0.8 0.88 0.8; 0.9 0.8 0.94]
    JBB_11 = SmoothLennardJones(ϵ[1,1], σ[1,1])
    JBB_12 = SmoothLennardJones(ϵ[1,2], σ[1,2])
    JBB_13 = SmoothLennardJones(ϵ[1,3], σ[1,3])
    JBB_22 = SmoothLennardJones(ϵ[2,2], σ[2,2])
    JBB_23 = SmoothLennardJones(ϵ[2,3], σ[2,3])
    JBB_33 = SmoothLennardJones(ϵ[3,3], σ[3,3])
    return SMatrix{3, 3, typeof(JBB_11), 9}([JBB_11 JBB_12 JBB_13; JBB_12 JBB_22 JBB_23; JBB_13 JBB_23 JBB_33])
    #return [JBB_11 JBB_12 JBB_13; JBB_12 JBB_22 JBB_23; JBB_13 JBB_23 JBB_33]
end

###############################################################################

struct GeneralKG{T<:AbstractFloat} <: DiscreteModel
    name::String
    ϵ::T
    ϵ4::T
    σ2::T
    shift::T
    k::T
    kr02::T
    r02::T
    rcut::T
    rcut2::T

end

function GeneralKG(ϵ, σ, k, r0; rcut=2^(1 / 6) * σ)
    name = "GeneralKG"
    r02 = r0 ^ 2
    rcut2 = rcut ^ 2
    kr02 = - k * r02 / 2
    σ2 = σ ^ 2
    shift = lennard_jones(rcut2, 4ϵ, σ2)
    return GeneralKG(name, ϵ, 4ϵ, σ2, shift, k, kr02, r02, rcut, rcut2)
end

@inline function potential(r2::T, model::GeneralKG) where T<:AbstractFloat
    return lennard_jones(r2, model.ϵ4, model.σ2) - model.shift
end

@inline function bond_potential(r2::T, model::GeneralKG) where T<:AbstractFloat
    return r2 ≤ model.r02 ? fene(r2, model.kr02, model.r02) : Inf
end

cutoff(model::DiscreteModel) = model.rcut
cutoff2(model::DiscreteModel) = model.rcut2

function Trimer()
    ϵ = SMatrix{3, 3, Float64}([1.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 1.0])
    σ = SMatrix{3, 3, Float64}([0.9 0.95 1.0; 0.95 1.0 1.05; 1.0 1.05 1.1])
    k = SMatrix{3, 3, Float64}([0.0 33.241 30.0; 33.241 0.0 27.210884; 30.0 27.210884 0.0])
    r0 = SMatrix{3, 3, Float64}([0.0 1.425 1.5; 1.425 0.0 1.575; 1.5 1.575 0.0])
    KG_11 = GeneralKG(ϵ[1,1], σ[1,1], k[1,1], r0[1,1])
    KG_12 = GeneralKG(ϵ[1,2], σ[1,2], k[1,2], r0[1,2])
    KG_13 = GeneralKG(ϵ[1,3], σ[1,3], k[1,3], r0[1,3])
    KG_22 = GeneralKG(ϵ[2,2], σ[2,2], k[2,2], r0[2,2])
    KG_23 = GeneralKG(ϵ[2,3], σ[2,3], k[2,3], r0[2,3])
    KG_33 = GeneralKG(ϵ[3,3], σ[3,3], k[3,3], r0[3,3])
    return SMatrix{3, 3, typeof(KG_11), 9}([KG_11 KG_12 KG_13; KG_12 KG_22 KG_23; KG_13 KG_23 KG_33])
end
###############################################################################