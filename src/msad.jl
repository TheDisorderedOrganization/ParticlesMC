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

    vals, vecs = eigen(R) # eigenvalues and eigenvectors of R
    idx = argmin(abs.(real.(vals) .- 1)) # find the idx corresponding to eigenvalue = 1
    n = real.(vecs[:, idx]) # extract the rotation axis vector
    n = SVector{3,T}(n[1], n[2], n[3])

    n_skew = SMatrix{3,3,T}(
         0,    n[3], -n[2],
        -n[3],  0,    n[1],
         n[2], -n[1],  0
    ) # skew matrix for n vector

    # Rodrigues formula
    sin_θ = clamp(-tr(n_skew * R) / 2, -one(T), one(T))

    # theta 
    θ     = atan(sin_θ, cos_θ)

    return θ * n
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
    states::Vector{MSADState{T}}
    theta_T::T
    compute_schedule::Vector{Int}
    output_schedule::Vector{Int}
    paths_integral::Vector{String}
    paths_thresh::Vector{String}
    files_integral::Vector{IOStream}
    files_thresh::Vector{IOStream}
end

function MSADTracker(chains;
                     theta_T::Float64=π/4,
                     output_schedule::Vector{Int}=Int[],
                     path::String=".",
                     scheduler::Vector{Int}=Int[],
                     kwargs...)

    compute_schedule = scheduler

    if !isempty(compute_schedule) && !isempty(output_schedule)
        @assert all(t in compute_schedule for t in output_schedule) """
        output_schedule contains steps not in compute_schedule.
        You cannot write output at a step where make_step! was not called.
        """
    end

    n      = length(chains)
    states = [MSADState{Float64}() for _ in 1:n]
    dirs   = [joinpath(path, "chains", "$c") for c in 1:n]

    paths_integral = [joinpath(d, "phi_integral.dat") for d in dirs]
    paths_thresh   = [joinpath(d, "phi_thresh.dat")   for d in dirs]

    files_integral = Vector{IOStream}(undef, n)
    files_thresh   = Vector{IOStream}(undef, n)

    return MSADTracker{Float64}(states, theta_T, compute_schedule,
                                output_schedule, paths_integral, paths_thresh,
                                files_integral, files_thresh)
end

##### Function to write phi frame #####
##########                   ##########

function write_phi_frame(file::IOStream, t::Int, N_mol::Int,
                         phis::Vector{<:SVector})
    println(file, N_mol)
    println(file, "t=$t")
    for m in 1:N_mol
        println(file, "$m $(phis[m][1]) $(phis[m][2]) $(phis[m][3])")
    end
    flush(file)
end

##### Initialisation of the simulation #####
# called once at t = 0, set the system, the mutable struct and the storing paths 
##########                        ##########

function Arianna.initialise(algorithm::MSADTracker, simulation::Simulation)

    for c in eachindex(simulation.chains)
        system = simulation.chains[c]
        state  = algorithm.states[c]
        T      = typeof(system.temperature)

        # create output directory if it doesn't exist
        mkpath(dirname(algorithm.paths_integral[c]))

        # open output files
        algorithm.files_integral[c] = open(algorithm.paths_integral[c], "w")
        algorithm.files_thresh[c]   = open(algorithm.paths_thresh[c],   "w")

        # compute initial body frames
        R_all = get_all_body_frames(system)
        N_mol = system.Nmol

        # fill state
        state.R0           = copy(R_all)
        state.R_prev       = copy(R_all)
        state.phi_integral = [zero(SVector{3,T}) for _ in 1:N_mol]
        state.R_ref_thresh = copy(R_all)
        state.phi_acc      = [zero(SVector{3,T}) for _ in 1:N_mol]
        state.initialized  = true

        # write t=0, all phi vectors are zero
        write_phi_frame(algorithm.files_integral[c], 0, N_mol,
                        state.phi_integral)
        write_phi_frame(algorithm.files_thresh[c],   0, N_mol,
                        [zero(SVector{3,T}) for _ in 1:N_mol])
    end
end

##### Finalise the simulation #####
# close output to not keep unwritten data in memory
##########              ##########

function Arianna.finalise(tracker::MSADTracker, ::Simulation)
    for c in eachindex(tracker.files_integral)
        close(tracker.files_integral[c])
        close(tracker.files_thresh[c])
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
                state.phi_acc[m]      = state.phi_total[m]
                state.R_ref_thresh[m] = R_all[m]
            end
        end

        state.R_prev .= R_all

        # output if time match output schedule
        if t in algorithm.output_schedule

            # Integral
            write_phi_frame(algorithm.files_integral[c],t,N_mol,state.phi_integral)

            # Threshold
            write_phi_frame(algorithm.files_integral[c],t,N_mol,state.phi_total)
        end
    end
end