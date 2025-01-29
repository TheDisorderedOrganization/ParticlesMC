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

cutoff(spi, spj, model::DiscreteModel) = model.rcut[spi, spj]
cutoff2(spi, spj, model::DiscreteModel) = model.rcut2[spi, spj]

"""
    abstract type ContinuousModel end

An abstract base type representing a general interaction model with a continuous size polidispersity.
"""
abstract type ContinuousModel <: Model end

cutoff(si, sj, model::ContinuousModel) = model.rcut * (si + sj) / 2

## POTENTIALS
harmonic_spheres(r2, epsilon, sigma2) = epsilon * (1 - r2 / sigma2)

inverse_power(r2, epsilon, sigma2, n::Int) = epsilon * (sigma2 / r2)^(n/2)

lennard_jones(r2, epsilon, sigma2) = 4 * epsilon * ((sigma2 / r2)^6 - (sigma2 / r2)^3)

wca(r2, epsilon, sigma2) = lennard_jones(r2, epsilon, sigma2) + epsilon

fene(r2, kr02, r02) = kr02 * log(1 - r2 / r02)

###############################################################################
"""
    struct SoftSpheres{T<:AbstractArray} <: DiscreteModel

A struct representing a soft-sphere interaction model where particles interact via an inverse power-law potential.

# Fields
- `name::String`: A descriptive name for the model.
- `epsilon::T`: Energy scale matrix for particle interactions between different types.
- `sigma::T`: Characteristic length matrix for particle interactions between different types.
- `n::Int`: Exponent controlling the softness of the interaction.
- `rcut::T`: Cutoff distance matrix beyond which interactions are neglected.
- `shift::T`: Shift matrix to ensure the potential is zero at the cutoff distance.
"""
struct SoftSpheres{T<:AbstractArray} <: DiscreteModel
    name::String
    epsilon::T
    sigma::T
    sigma2::T
    n::Int
    rcut::T
    rcut2::T
    shift::T
end

function SoftSpheres(epsilon, sigma, n; rcut=2.5*sigma, name="SoftSpheres$(sigma[2,2]/sigma[1,1])")
    shift = [inverse_power(rcut[spi, spj], epsilon[spi, spj], sigma[spi, spj], n) for spi in 1:2, spj in 1:2]
    sigma2 = sigma .^ 2
    rcut2 = rcut .^ 2
    return SoftSpheres(name, epsilon, sigma, sigma2, n, rcut, rcut2, shift)
end

function BHHP()
    epsilon = [1.0 1.0; 1.0 1.0]
    sigma = [1.0 1.2; 1.2 1.4]
    return SoftSpheres(epsilon, sigma, 12; rcut=2.5*sigma, name="BHHP")
end

function potential(r2, spi, spj, model::SoftSpheres)
    return inverse_power(r2, model.epsilon[spi, spj], model.sigma2[spi, spj], model.n) - model.shift[spi, spj]
end

###############################################################################
"""
    struct KobAndersen{T<:AbstractArray} <: DiscreteModel

A struct representing the KobAndersen interaction model where particles interact via a Lennard-Jones potential.

# Fields
- `name::String`: A descriptive name for the model.
- `epsilon::T`: Energy scale matrix for particle interactions between different types.
- `sigma::T`: Characteristic length matrix for particle interactions between different types.
- `rcut::T`: Cutoff distance matrix beyond which interactions are neglected.
- `shift::T`: Shift matrix to ensure the potential is zero at the cutoff distance.
"""
struct KobAndersen{T<:AbstractArray} <: DiscreteModel
    name::String
    epsilon::T
    sigma::T
    sigma2::T
    rcut::T
    rcut2::T
    shift::T
end

function KobAndersen(epsilon, sigma; rcut=2.5*sigma, name="KA$(size(sigma)[1]-1)")
    shift = [lennard_jones(rcut[spi, spj], epsilon[spi, spj], sigma[spi, spj]) for spi in axes(sigma, 1), spj in axes(sigma, 1)]
    sigma2 = sigma .^ 2
    rcut2 = rcut .^ 2
    return KobAndersen(name, epsilon, sigma, sigma2, rcut, rcut2, shift)
end

function KobAndersen()
    epsilon = [1.0 1.5; 1.5 0.5]
    sigma = [1.0 0.8; 0.8 0.88]
    return KobAndersen(epsilon, sigma; rcut=2.5*sigma, name="KobAndersen")
end

function potential(r2, spi, spj, model::KobAndersen)
    return lennard_jones(r2, model.epsilon[spi, spj], model.sigma2[spi, spj]) - model.shift[spi, spj]
end

###############################################################################

struct JBB{V<:AbstractArray, T<:AbstractFloat} <: DiscreteModel
    name::String
    epsilon::V
    sigma::V
    sigma2::V
    sigma4::V
    rcut::V
    rcut2::V
    C0::T
    C2::T
    C4::T
end

function JBB()
    name = "JBB"
    epsilon = SMatrix{3,3,Float64}([1.0 1.5 0.75; 1.5 0.5 1.5; 0.75 1.5 0.75])
    sigma = SMatrix{3,3,Float64}([1.0 0.8 0.9; 0.8 0.88 0.8; 0.9 0.8 0.94])
    rcut = 2.5 * sigma
    C0 = 0.04049023795
    C2 = -0.00970155098
    C4 = 0.00062012616
    sigma2 = sigma .^ 2
    sigma4 = sigma2 .^ 2
    rcut2 = rcut .^ 2
    return JBB(name, epsilon, sigma, sigma2, sigma4, rcut, rcut2, C0, C2, C4)
end

function potential(r2, spi, spj, model::JBB)
    epsilon, sigma2, sigma4 = model.epsilon[spi, spj], model.sigma2[spi, spj], model.sigma4[spi, spj]
    lj = lennard_jones(r2, epsilon, sigma2)
    shift = model.C0 + model.C2 * r2 / sigma2 + model.C4 * r2^2 / sigma4
    return lj + 4 * epsilon * shift
end

###############################################################################

struct GeneralKG{V<:AbstractArray} <: DiscreteModel
    name::String
    epsilon::V
    sigma2::V
    k::V
    r02::V
    kr02::V
    rcut::V
    rcut2::V
end

function GeneralKG(epsilon, sigma, k, r0; rcut=2^(1 / 6) * sigma)
    name = "GeneralKG"
    sigma2 = sigma .^ 2
    r02 = r0 .^ 2
    rcut2 = rcut .^ 2
    kr02 = -0.5 .* k .* r02 
    return GeneralKG(name, epsilon, sigma2, k, r02, kr02, rcut, rcut2)
end

function potential(r, spi, spj, model::GeneralKG)
    return wca(r, model.epsilon[spi, spj], model.sigma2[spi, spj])
end

function bond_potential(r2, spi, spj, model::GeneralKG)
    r02 = model.r02[spi, spj]
    if r2 ≤ r02
        res = fene(r2, model.kr02[spi, spj], r02)
    else
        res = Inf
    end
    return res
end

cutoff(spi, spj, model::GeneralKG) = model.rcut[spi, spj]
cutoff2(spi, spj, model::GeneralKG) = @inbounds model.rcut2[spi, spj]

###############################################################################

nothing
