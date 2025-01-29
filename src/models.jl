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

"""
    abstract type ContinuousModel end

An abstract base type representing a general interaction model with a continuous size polidispersity.
"""
abstract type ContinuousModel <: Model end

cutoff(si, sj, model::ContinuousModel) = model.rcut * (si + sj) / 2

## POTENTIALS
harmonic_spheres(r, epsilon, sigma) = epsilon * (1 - r / sigma)^2

inverse_power(r, epsilon, sigma, n::Int) = epsilon * (sigma / r)^n

lennard_jones(r, epsilon, sigma) = 4 * epsilon * ((sigma / r)^12 - (sigma / r)^6)

wca(r, epsilon, sigma) = lennard_jones(r, epsilon, sigma) + epsilon

fene(r, k, r0) = -0.5 * k * r0^2 * log(1 - (r / r0)^2)

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
    n::Int
    rcut::T
    shift::T
end

function SoftSpheres(epsilon, sigma, n; rcut=2.5*sigma, name="SoftSpheres$(sigma[2,2]/sigma[1,1])")
    shift = [inverse_power(rcut[spi, spj], epsilon[spi, spj], sigma[spi, spj], n) for spi in 1:2, spj in 1:2]
    return SoftSpheres(name, epsilon, sigma, n, rcut, shift)
end

function BHHP()
    epsilon = [1.0 1.0; 1.0 1.0]
    sigma = [1.0 1.2; 1.2 1.4]
    return SoftSpheres(epsilon, sigma, 12; rcut=2.5*sigma, name="BHHP")
end

function potential(r, spi, spj, model::SoftSpheres)
    return inverse_power(r, model.epsilon[spi, spj], model.sigma[spi, spj], model.n) - model.shift[spi, spj]
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
    rcut::T
    shift::T
end

function KobAndersen(epsilon, sigma; rcut=2.5*sigma, name="KA$(size(sigma)[1]-1)")
    shift = [lennard_jones(rcut[spi, spj], epsilon[spi, spj], sigma[spi, spj]) for spi in axes(sigma, 1), spj in axes(sigma, 1)]
    return KobAndersen(name, epsilon, sigma, rcut, shift)
end

function KobAndersen()
    epsilon = [1.0 1.5; 1.5 0.5]
    sigma = [1.0 0.8; 0.8 0.88]
    return KobAndersen(epsilon, sigma; rcut=2.5*sigma, name="KobAndersen")
end

function potential(r, spi, spj, model::KobAndersen)
    return lennard_jones(r, model.epsilon[spi, spj], model.sigma[spi, spj]) - model.shift[spi, spj]
end

###############################################################################

struct JBB{V<:AbstractArray, T<:AbstractFloat} <: DiscreteModel
    name::String
    epsilon::V
    sigma::V
    rcut::V
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
    return JBB(name, epsilon, sigma, rcut, C0, C2, C4)
end

function potential(r, spi, spj, model::JBB)
    epsilon_ij, sigma_ij = model.epsilon[spi, spj], model.sigma[spi, spj]
    lj = lennard_jones(r, epsilon_ij, sigma_ij)
    shift = model.C0 + model.C2 * (r / sigma_ij)^2 + model.C4 * (r / sigma_ij)^4
    return lj + 4 * epsilon_ij * shift
end

###############################################################################

struct GeneralKG{V<:AbstractArray} <: DiscreteModel
    name::String
    epsilon::V
    sigma::V
    k::V
    r0::V
    rcut::V
end

function GeneralKG(epsilon, sigma, k, r0; rcut=2^(1 / 6) * sigma)
    name = "GeneralKG"
    return GeneralKG(name, epsilon, sigma, k, r0, rcut)
end

function potential(r, spi, spj, model::GeneralKG)
    return wca(r, model.epsilon[spi[1], spj[1]], model.sigma[spi[1], spj[1]])
end

function bond_potential(r, spi, spj, model::GeneralKG)
    if r ≤ model.r0[spi[1], spj[1]]
        res = fene(r, model.k[spi[1], spj[1]], model.r0[spi[1], spj[1]])
    else
        res = Inf
    end
    return res
end

cutoff(spi, spj, model::GeneralKG) = model.rcut[spi[1], spj[1]]

###############################################################################

nothing
