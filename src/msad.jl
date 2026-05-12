# msad.jl
#
# On-the-fly Mean Squared Angular Displacement (MSAD) tracker.
# Implements three methods:
#   - Integral  : accumulates incremental rotation vectors every step
#   - Threshold : relays the reference frame when rotation exceeds θ_T
#   - Euler     : compares directly to initial frame, computed on log schedule

using LinearAlgebra
using StaticArrays

##### Basic operations for rotation inspired fromo our python post processing #####
##########                                                               ##########

function body_frame(r1::SVector{3,T}, r2::SVector{3,T}, r3::SVector{3,T}, L::T) where {T}
    e1 = r2 - r1
    v  = r3 - r1
    # minimum image convention
    e1 = e1 - round.(e1 / L) .* L
    v  =  v - round.( v / L) .* L
    # Gram-Schmidt basis right handed
    e1 = e1 / norm(e1)
    e2 = cross(e1, v)
    e2 = e2 / norm(e2)
    e3 = cross(e1, e2)
    # hcat stacks column vectors → 3×3 matrix
    return hcat(e1, e2, e3)   # SMatrix{3,3,T}
end

function get_all_body_frames(system::Molecules)
    L   = system.box[1]                                    # cubic box
    T   = typeof(L)
    R   = Vector{SMatrix{3,3,T,9}}(undef, system.Nmol)    # pre-allocate
    pos = system.position
    @inbounds for m in 1:system.Nmol                    # iterate over the number of molecules
        s    = system.start_mol[m]                         # first bead index for the considered molecule
        R[m] = body_frame(pos[s], pos[s+1], pos[s+2], L)
    end
    return R
end

function rotation_vector(R::SMatrix{3,3,T}) where {T}
    cos_θ = clamp((tr(R) - 1) / 2, -one(T), one(T))
    # skew-symmetric part encodes sin(θ)*n
    skew  = (R - R') / 2
    sin_n = SVector(skew[3,2], skew[1,3], skew[2,1]) # julia starts index at 1
    sin_θ = norm(sin_n)
    θ     = atan(sin_θ, cos_θ)
    return sin_θ > 1e-10 ? (θ / sin_θ) * sin_n : zero(SVector{3,T}) ## might need to add the theta = pi case (unprobable for MC steps)
end

##### Definition of the MSAD State ie: the accumulating variables #####
# mutable because we need to update it 
##########                                                   ##########

mutable struct MSADState{T}
    R0::Vector{SMatrix{3,3,T,9}}            # (Nmol,3,3) initial body frames at t=0
    R_prev::Vector{SMatrix{3,3,T,9}}        # (Nmol,3,3) body frames at previous time 
    phi_integral::Vector{SVector{3,T}}      # (Nmol,3) Accumulated vector for integral method
    R_ref_thresh::Vector{SMatrix{3,3,T,9}}  # (Nmol,3,3) reference body frames for each molecule in the threshold method
    phi_acc::Vector{SVector{3,T}}           # (Nmol,3) Accumulated vector for threshold method
    initialized::Bool
end

# Start empty, filled during initialise when we first see the system
MSADState{T}() where {T} = MSADState{T}([], [], [], [], [], false)

##### The MSAD computation algorithm #####
##########                      ##########

mutable struct MSADTracker{T} <: AriannaAlgorithm
    states::Vector{MSADState{T}}    # the states of the considered chain
    theta_T::T                      # threshold angle in radians
    output_schedule::Vector{Int}    # schedule for writing MSAD to disk 
    paths::Vector{String}
    files::Vector{IOStream}
end

function MSADTracker(chains, theta_T::T, output_schedule::Vector{Int}, path::String) where {T}
    n = length(chains)
    states = [MSADState{T}() for _ in 1:n] # one empty state per chain filled during initialise
    
    dirs  = [joinpath(path, "chains", "$c") for c in 1:n]
    paths = [joinpath(d, "msad.dat") for d in dirs]
    files = Vector{IOStream}(undef, n)
    return MSADTracker{T}(states, theta_T, output_schedule, paths, files)
end