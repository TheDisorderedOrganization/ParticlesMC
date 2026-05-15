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
    compute_schedule::Vector{Int}   # schedule for computing rotation for integral and threshold methods
    output_schedule::Vector{Int}    # schedule for writing MSAD to disk 
    paths::Vector{String}
    files::Vector{IOStream}
end

function MSADTracker(chains, theta_T::T, compute_schedule::Vector{Int}, output_schedule::Vector{Int}, path::String) where {T}

    @assert all(t in compute_schedule for t in output_schedule) # output_schedule must be a subset of compute_schedule

    n = length(chains)
    states = [MSADState{T}() for _ in 1:n] # one empty state per chain filled during initialise

    dirs  = [joinpath(path, "chains", "$c") for c in 1:n]
    paths = [joinpath(d, "msad.dat") for d in dirs]
    files = Vector{IOStream}(undef, n)
    return MSADTracker{T}(states, theta_T, compute_schedule, output_schedule, paths, files)
end

##### Initialisation of the simulation #####
# called once at t = 0, set the system, the mutable struct and the storing paths 
##########                        ##########

function Arianna.initialise(algorithm::MSADTracker, simulation::Simulation)

    for c in eachindex(simulation.chains)
        system = simulation.chains[c]
        state  = algorithm.states[c]
        T      = typeof(system.temperature)   # get the float type from system so one can choose Float64 or 32 (generic)

        # create output directory if it doesn't exist
        mkpath(dirname(algorithm.paths[c]))

        # open output file — write header
        algorithm.files[c] = open(algorithm.paths[c], "w")
        println(algorithm.files[c], "# t  msad_euler  msad_integral  msad_thresh")

        # compute initial body frames
        R_all = get_all_body_frames(system)
        N_mol = system.Nmol

        # fill state, copy. so they are different pointers!
        state.R0           = copy(R_all)
        state.R_prev       = copy(R_all)
        state.phi_integral = [zero(SVector{3,T}) for _ in 1:N_mol]
        state.R_ref_thresh = copy(R_all)
        state.phi_acc      = [zero(SVector{3,T}) for _ in 1:N_mol]
        state.initialized  = true

        # write t=0, all MSAD values are zero at start
        println(algorithm.files[c], "0 0.0 0.0 0.0")
        flush(algorithm.files[c])
    end

end

##### Finalise th simulation #####
# close output to not keep unwritten data in memory
##########              ##########

function Arianna.finalise(algorithm::MSADTracker, ::Simulation)
    for f in algorithm.files
        close(f)
    end
end

##### Make a step in simulation #####
##########                 ##########

function Arianna.make_step!(simulation::Simulation, algorithm::MSADTracker)
    t = simulation.t

    for c in eachindex(simulation.chains)
        system = simulation.chains[c]
        state  = algorithm.states[c]
        N_mol  = system.Nmol

        # compute current body frames shared by all three methods
        R_all = get_all_body_frames(system)

        # Integral update
        for m in 1:N_mol
            dR = state.R_prev[m]' * R_all[m]
            state.phi_integral[m] = state.phi_integral[m] + rotation_vector(dR)
        end

        # Threshold update
        phi_current = Vector{SVector{3,eltype(algorithm.theta_T)}}(undef, N_mol)
        phi_total   = Vector{SVector{3,eltype(algorithm.theta_T)}}(undef, N_mol)
        for m in 1:N_mol
            dR             = state.R_ref_thresh[m]' * R_all[m]
            phi_current[m] = rotation_vector(dR)
            phi_total[m]   = state.phi_acc[m] + phi_current[m]
            if norm(phi_current[m]) >= algorithm.theta_T
                state.phi_acc[m]      = phi_total[m]
                state.R_ref_thresh[m] = R_all[m]
            end
        end

        state.R_prev .= R_all

        # output is time match output schedule
        if t in algorithm.output_schedule
            # Euler only computed here not accumulated
            msad_euler = sum(
                norm(rotation_vector(state.R0[m]' * R_all[m]))^2
                for m in 1:N_mol
            ) / N_mol

            # Integral
            msad_integral = sum(
                norm(state.phi_integral[m])^2
                for m in 1:N_mol
            ) / N_mol

            # Threshold
            msad_thresh = sum(
                norm(phi_total[m])^2
                for m in 1:N_mol
            ) / N_mol

            println(algorithm.files[c], "$t $msad_euler $msad_integral $msad_thresh")
            flush(algorithm.files[c])
        end
    end
end