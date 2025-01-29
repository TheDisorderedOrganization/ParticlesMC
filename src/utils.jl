
using ConcreteStructs

function MonteCarlo.delta_log_target_density(e1, e2, system::Particles)
    return -(e2 - e1) ./ system.temperature
end

function nearest_image_distance(xi::T, xj::T, L::T) where {T<:AbstractFloat}
    dx = xi - xj
    return dx - round(dx / L) * L
end

function nearest_image_distance(xi::T, xj::T, box::T) where {T<:AbstractArray}
    return map((xi, xj, L) -> nearest_image_distance(xi, xj, L), xi, xj, box)
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

nothing