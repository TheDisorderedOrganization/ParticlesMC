# rotation.jl
#
# Compute on-the-fly rotation tracker.
# Two algorithms :
# 1. ComputeRotation : updates system.Φ (simulation.observable)
# 2. StorePhiTrajectory : writes system.Φ to disk
# system.Φ[k][m] = rotation vector for molecule m under θ_T[k]

using LinearAlgebra
using StaticArrays

##### Basic operations for rotation #####
##########                          ##########

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

function get_all_body_frames!(R_bodyframe, system::Molecules)
    L   = system.box[1]
    pos = system.position
    @inbounds for m in 1:system.Nmol
        s = system.start_mol[m]
        R_bodyframe[m] = body_frame(pos[s], pos[s+1], pos[s+2], L)
    end
end

function rotation_vector(R::SMatrix{3,3,T}) where {T}
    cos_θ = clamp((tr(R) - 1) / 2, -one(T), one(T))

    ax = R[3,2] - R[2,3]
    ay = R[1,3] - R[3,1]
    az = R[2,1] - R[1,2]
    sin_θ = sqrt(ax^2 + ay^2 + az^2) / 2

    if sin_θ > sqrt(eps(T))              # guard sin_θ ≈ 0
        θ      = atan(sin_θ, cos_θ)
        inv_2s = θ / (2 * sin_θ)
        return SVector{3,T}(ax * inv_2s, ay * inv_2s, az * inv_2s)
    
    elseif cos_θ > 0                     # guard θ ≈ 0, near identity
        return zero(SVector{3,T})

    else                                  # guard θ ≈ π
        if R[1,1] >= R[2,2] && R[1,1] >= R[3,3]
            nx = sqrt(max((R[1,1] + 1) / 2, zero(T)))
            n  = SVector{3,T}(nx, (R[1,2] + R[2,1]) / (4nx), (R[1,3] + R[3,1]) / (4nx))
        elseif R[2,2] >= R[3,3]
            ny = sqrt(max((R[2,2] + 1) / 2, zero(T)))
            n  = SVector{3,T}((R[1,2] + R[2,1]) / (4ny), ny, (R[2,3] + R[3,2]) / (4ny))
        else
            nz = sqrt(max((R[3,3] + 1) / 2, zero(T)))
            n  = SVector{3,T}((R[1,3] + R[3,1]) / (4nz), (R[2,3] + R[3,2]) / (4nz), nz)
        end
        θ = atan(sin_θ, cos_θ)
        return θ * n / norm(n)
    end
end

##### Definition of the RotationState ie: the accumulating variables #####
# mutable because we need to update it 
##########                                                   ##########

mutable struct RotationState{T}
    R_ref::Vector{Vector{SMatrix{3,3,T,9}}} 
    R_bodyframe::Vector{SMatrix{3,3,T,9}} 
    Φ_acc::Vector{Vector{SVector{3,T}}}       
    initialized::Bool
end

# Start empty, filled during initialise when we first see the system
RotationState{T}() where {T} = RotationState{T}([], [], [], false)

##### The Rotation computation algorithm #####
##########                      ##########

mutable struct ComputeRotation{T} <: AriannaAlgorithm
    θ_T::Vector{T}                  # one per method
    states::Vector{RotationState{T}}    # one per chain
end

function ComputeRotation(chains;
                     θ_T::Vector{Float64}=[π/4],
                     scheduler::Vector{Int}=Int[],
                     kwargs...)
    n = length(chains) # number of chains
    states = [RotationState{Float64}() for _ in 1:n] # n rotation states independent
    return ComputeRotation{Float64}(θ_T, states)
end
    
##### Initialisation of the simulation #####
# called once at t = 0, set the system, the mutable struct and the storing paths 
##########                        ##########

function Arianna.initialise(algorithm::ComputeRotation, simulation::Simulation)
    for c in eachindex(simulation.chains)
        system = simulation.chains[c]
        state  = algorithm.states[c]
        T      = typeof(system.temperature)
        N_mol  = system.Nmol
        n_θ    = length(algorithm.θ_T)

        state.R_bodyframe = Vector{SMatrix{3,3,T,9}}(undef, N_mol)  # allocate once
        get_all_body_frames!(state.R_bodyframe, system)              # fill state.R_bodyframe

        if simulation.t_start > 0
            _, _, _, _, R_ref, Φ_acc, Φ = read_phi_frame(joinpath(simulation.path,"chains", "$c", "lastphiframe.dat"))
            state.R_ref = R_ref
            state.Φ_acc = Φ_acc
            resize!(system.Φ, n_θ)
            for k in 1:n_θ
                system.Φ[k] = Φ[k]
            end
        
        else 
            state.R_ref       = [copy(state.R_bodyframe) for _ in 1:n_θ]
            state.Φ_acc       = [[zero(SVector{3,T}) for _ in 1:N_mol] for _ in 1:n_θ]
            resize!(system.Φ, n_θ)
            for k in 1:n_θ
                system.Φ[k] = [zero(SVector{3,T}) for _ in 1:N_mol]
            end
        end
        state.initialized = true
    end
end

##### Make a step in simulation #####
##########                 ##########

function Arianna.make_step!(simulation::Simulation, algorithm::ComputeRotation)
    tcollect(eachindex(simulation.chains) |> Transducers.Map(c -> begin
        system = simulation.chains[c]
        state  = algorithm.states[c]
        N_mol  = system.Nmol
        get_all_body_frames!(state.R_bodyframe, system)   # no new vector created updates state.R_bodyframe
        for (k, θ_T) in enumerate(algorithm.θ_T)
            for m in 1:N_mol
                dR        = state.R_ref[k][m]' * state.R_bodyframe[m]
                Φ_current = rotation_vector(dR)
                Φ_total   = state.Φ_acc[k][m] + Φ_current
                if norm(Φ_current) > θ_T
                    state.Φ_acc[k][m] = Φ_total
                    state.R_ref[k][m] = state.R_bodyframe[m]
                end
                system.Φ[k][m] = Φ_total
            end
        end
    end))
end

function Arianna.finalise(::ComputeRotation, ::Simulation) end

##### Read a checkpoint ɸ frame #####
##########                 ##########

function read_phi_frame(path)
    open(path) do file
        t = parse(Int, split(readline(file), "=")[2])
        N_mol = parse(Int, split(readline(file), "=")[2])
        n_θ = parse(Int, split(readline(file), "=")[2])
        θ_T = parse.(Float64, split(split(readline(file), "=")[2], ","))

        R_ref = Vector{Vector{SMatrix{3,3,Float64,9}}}(undef, n_θ)
        Φ_acc = Vector{Vector{SVector{3,Float64}}}(undef, n_θ)
        Φ     = Vector{Vector{SVector{3,Float64}}}(undef, n_θ)

        for k in 1:n_θ
            R_ref[k] = read_block(file, N_mol, row_to_matrix)
            Φ_acc[k] = read_block(file, N_mol, row_to_vector)
            Φ[k] = read_block(file, N_mol, row_to_vector)
        end
    return (t, N_mol, n_θ, θ_T, R_ref, Φ_acc, Φ)
    end
end

function read_block(file, N_mol, builder)
    readline(file) # skip first # line
    return [builder(split(readline(file), " ")) for _ in 1:N_mol]
end

function row_to_matrix(parts)
    return SMatrix{3,3,Float64}(parse.(Float64, parts[2:10]))
end

function row_to_vector(parts)
    return SVector{3,Float64}(parse.(Float64, parts[2:4]))
end

##### Store rotation vector trajectory #####
##########                        ##########

function write_phi_frame(file::IOStream, t::Int, N_mol::Int, Φs::Vector{<:SVector})
    println(file, N_mol)
    println(file, "t=$t")
    for m in 1:N_mol
        println(file, "$m $(Φs[m][1]) $(Φs[m][2]) $(Φs[m][3])")
    end
    #flush(file) not flushing at everyframe
end

struct StorePhiTrajectories <: AriannaAlgorithm
    dirs::Vector{String}                # base dirs per chain
    paths::Vector{Vector{String}}       # paths[c][k]
    files::Vector{Vector{IOStream}}     # files[c][k]
    store_first::Bool
    store_last::Bool

    function StorePhiTrajectories(chains, path; store_first::Bool=true, store_last::Bool=false)
        dirs  = joinpath.(path, "chains", ["$(c)" for c in eachindex(chains)])
        mkpath.(dirs)
        n     = length(chains)
        paths = [String[] for _ in 1:n]
        files = [IOStream[] for _ in 1:n]
        return new(dirs, paths, files, store_first, store_last)
    end
end

function StorePhiTrajectories(chains; path=missing, store_first=true, store_last=false, kwargs...)
    return StorePhiTrajectories(chains, path; store_first=store_first, store_last=store_last)
end

function Arianna.initialise(algorithm::StorePhiTrajectories, simulation::Simulation)
    simulation.verbose && println("Opening Φ trajectory files...")
    for c in eachindex(simulation.chains)
        system = simulation.chains[c]
        writing_mode = simulation.t_start > 0 ? "a" : "w"
        n_θ    = length(system.Φ)
        algorithm.paths[c] = [joinpath(algorithm.dirs[c], "phitrajectories_$k.dat")
                               for k in 1:n_θ]
        algorithm.files[c] = open.(algorithm.paths[c], writing_mode)
    end
    (simulation.t_start == 0 && algorithm.store_first) && Arianna.make_step!(simulation, algorithm)
end

function Arianna.make_step!(simulation::Simulation, algorithm::StorePhiTrajectories)
    for c in eachindex(simulation.chains)
        system = simulation.chains[c]
        n_θ    = length(system.Φ)
        N_mol  = system.Nmol
        for k in 1:n_θ
            write_phi_frame(algorithm.files[c][k], simulation.t, N_mol, system.Φ[k])
        end
    end
end

function Arianna.finalise(algorithm::StorePhiTrajectories, simulation::Simulation)
    algorithm.store_last && make_step!(simulation, algorithm)
    simulation.verbose && println("Closing Φ trajectory files...")
    for files_c in algorithm.files
        close.(files_c)
    end
end

##### Store last Φ frame for checkpoint/restart #####
##########                                    ##########



struct StoreLastPhiFrame <: AriannaAlgorithm
    paths::Vector{String}   # one per chain: chains/c/rotation_checkpoint.dat

    function StoreLastPhiFrame(chains, path; kwargs...)
        dirs  = joinpath.(path, "chains", ["$(c)" for c in eachindex(chains)])
        mkpath.(dirs)
        paths = joinpath.(dirs, "lastphiframe.dat")
        return new(paths)
    end
end

function StoreLastPhiFrame(chains; path=missing, kwargs...)
    return StoreLastPhiFrame(chains, path)
end

function Arianna.finalise(algorithm::StoreLastPhiFrame, simulation::Simulation)
    # find ComputeRotation in simulation algorithms to access internal state
    compute_rot = nothing
    for algo in simulation.algorithms
        if algo isa ComputeRotation
            compute_rot = algo
            break
        end
    end
    @assert compute_rot !== nothing "StoreLastPhiFrame requires ComputeRotation in algorithm list"

    for c in eachindex(simulation.chains)
        system = simulation.chains[c]
        state  = compute_rot.states[c]
        N_mol  = system.Nmol
        n_θ    = length(compute_rot.θ_T)

        open(algorithm.paths[c], "w") do file
            # header
            println(file, "t=$(simulation.t)")
            println(file, "N_mol=$N_mol")
            println(file, "n_θ=$n_θ")
            println(file, "θ_T=$(join(compute_rot.θ_T, ','))")

            # R_ref and Φ_acc
            for k in 1:n_θ
                println(file, "# R_ref k=$k")
                for m in 1:N_mol
                    r = state.R_ref[k][m]
                    println(file, "$m $(r[1,1]) $(r[2,1]) $(r[3,1]) $(r[1,2]) $(r[2,2]) $(r[3,2]) $(r[1,3]) $(r[2,3]) $(r[3,3])")
                end
                println(file, "# Φ_acc k=$k")
                for m in 1:N_mol
                    p = state.Φ_acc[k][m]
                    println(file, "$m $(p[1]) $(p[2]) $(p[3])")
                end
                println(file, "# Φ k=$k")
                for m in 1:N_mol
                    p = system.Φ[k][m]
                    println(file, "$m $(p[1]) $(p[2]) $(p[3])")
                end
            end
        end
    end
    return nothing
end
