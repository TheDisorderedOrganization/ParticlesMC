
using ConcreteStructs, Distributions, Statistics, StaticArrays

function Arianna.unnormalised_log_target_density(e, system::Particles)
    return  - e / system.temperature
end

function Arianna.delta_log_target_density(e1, e2, system::Particles)
    return -(e2 - e1) ./ system.temperature
end

fold_back(x, box) = x .- fld.(x, box) .* box


function vector_1D(c1, c2, side_length)
    dx = c1 - c2
    return dx - round(dx / side_length) * side_length
end

@inline function vector(c1::SVector{N,T}, c2::SVector{N,T}, box::SVector{N,T}) where {N, T<:AbstractFloat}
    @inbounds return SVector{N,T}(ntuple(i -> vector_1D(c1[i], c2[i], box[i]), Val(N)))
end

@inline function nearest_image_distance_squared(xi::SVector{N,T}, xj::SVector{N,T},
                                                box::SVector{N,T}) where {N, T<:AbstractFloat}
    @inbounds dx = vector(xi, xj, box)
    return sum(abs2, dx)
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

Arianna.@callback function energy(system)
    return system.energy[1]
end

Arianna.@callback function chain_correlation(system)
    return compute_chain_correlation(system)
end

function volume_sphere(r, d::Int)
    d == 0 && return 1
    d == 1 && return 2r
    return 2π * r^2 * volume_sphere(r, d - 2) / d 
end