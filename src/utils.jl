
using ConcreteStructs, Distributions, Statistics, StaticArrays

function Arianna.delta_log_target_density(e1, e2, system::Particles)
    return -(e2 - e1) ./ system.temperature
end

fold_back(x, box) = x .- fld.(x, box) .* box


function vector_1D(c1, c2, side_length)
    dx = c1 - c2
    return dx - round(dx / side_length) * side_length
end

#Taken and adapted from Molly.jl
function vector(c1::SVector{2,T}, c2::SVector{2,T}, box::SVector{2,T}) where {T<:AbstractFloat}
    return @inbounds SVector(
        vector_1D(c1[1], c2[1], box[1]),
        vector_1D(c1[2], c2[2], box[2]),
    )
end


#Taken and adapted from Molly.jl
function vector(c1::SVector{3,T}, c2::SVector{3,T}, box::SVector{3,T}) where {T<:AbstractFloat}
    return @inbounds SVector(
        vector_1D(c1[1], c2[1], box[1]),
        vector_1D(c1[2], c2[2], box[2]),
        vector_1D(c1[3], c2[3], box[3]),
    )
end


function nearest_image_distance_squared(xi::SVector{N,T}, xj::SVector{N,T},
                                        box::SVector{N,T}) where {N,T<:AbstractFloat}
    dx = vector(xi, xj, box)
    return sum(abs2, dx)             # Sum of squared elements
end

struct SpeciesList
    sp_ids::Vector{Vector{Int}}
    sp_heads::Vector{Int}
end

function SpeciesList(species)
    available_species = sort(unique(species))
    ids = Vector{Vector{Int}}(undef, length(available_species))
    heads = zeros(Int, length(species))
    for k in eachindex(available_species)
        ids[k] = findall(isequal(available_species[k]), species)
        for i in eachindex(species)
            if species[i] == available_species[k]
                heads[i] = findfirst(isequal(i), ids[k])
            end
        end
    end
    return SpeciesList(ids, heads)
end

function callback_energy(simulation)
    return mean(system.energy[1] / length(system) for system in simulation.chains)
end


function volume_sphere(r, d::Int)
    d == 0 && return 1
    d == 1 && return 2r
    return 2π * r^2 * volume_sphere(r, d - 2) / d 
end