struct Atoms{D, V<:AbstractVector, M<:Model, C<:CellList, H, T<:AbstractFloat} <: Particles
    position::Vector{SVector{D,T}}
    species::V
    density::T
    temperature::T
    model::M
    d::Int
    N::Int
    box::SVector{D,T}
    local_energy::Vector{T}
    cell_list::C
    species_list::H
    cache::Vector{Tuple{Int,T}} # We only cache local energies for now
    particle_ids::Base.OneTo{Int}
end

function System(position, species, density::T, temperature::T, model::Model; list_type=EmptyList) where{T<:AbstractFloat}
    @assert length(position) == length(species)
    particle_ids = eachindex(position)
    N = length(particle_ids)
    d = length(Array(position)[1])
    box = @SVector fill(T((N / density)^(1 / d)), d)
    e_locals = zeros(T, N)
    maxcut = maximum([cutoff(spi, spj, model) for spi in species for spj in species])
    cell_list = list_type(box, maxcut, N)
    species_list = isa(species[1], Integer) ? SpeciesList(species) : nothing
    cache = Vector{Tuple{Int,T}}()
    system = Atoms(position, species, density, temperature, model, d, N, box, e_locals, cell_list, species_list, cache, particle_ids)
    build_cell_list!(system, system.cell_list)
    system.local_energy .= [local_energy(system, i, system.cell_list) for i in particle_ids]
    return system
end

function calc_energy_ij!(system::Atoms, i, j,  position_i, species_i)
    if i != j
        position_j, species_j = system.position[j], system.species[j]
        r = norm(nearest_image_distance(position_i, position_j, system.box))
        if r ≤ cutoff(species_i, species_j, system.model)
            energy_ij = potential(r, species_i, species_j, system.model)
            return energy_ij
        end
    end
    return zero(typeof(system.density))
end

function local_energy(system::Atoms, i, ::EmptyList)
    energy = zero(typeof(system.density))
    position_i, species_i = system.position[i], system.species[i]
    for j in system.particle_ids
        energy += calc_energy_ij!(system, i, j, position_i, species_i)
    end
    return energy
end

function local_energy(system::Atoms, i, cell_list::LinkedList)
    energy = zero(typeof(system.density))
    position_i, species_i = system.position[i], system.species[i]
    mc = get_cell(position_i, cell_list.cell)
    # Scan the neighbourhood of cell mc (including itself)
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(mc2, cell_list.ncells)
        # Scan atoms in cell c2
        j = cell_list.head[c2]
        while (j != -1)
            energy += calc_energy_ij!(system, i, j, position_i, species_i)
            j = cell_list.list[j]
        end
    end
    return energy
end

function destroy_particle!(system::Atoms, i, ::EmptyList)
    energy = zero(typeof(system.density))
    position_i, species_i = system.position[i], system.species[i]
    push!(system.cache, (i, system.local_energy[i]))
    # Loop over particles
    @inbounds for j in system.particle_ids
        push!(system.cache, (j, system.local_energy[j]))
        energy_ij = calc_energy_ij!(system, i, j, position_i, species_i)
        system.local_energy[j] -= energy_ij
        energy += energy_ij
    end
    # Update particle i local energy cache
    system.local_energy[i] = energy
    return energy
end

function create_particle!(system::Atoms, i, ::EmptyList)
    energy = zero(typeof(system.density))
    position_i, species_i = system.position[i], system.species[i]
    # Loop over particles
    @inbounds for j in system.particle_ids
        energy_ij = calc_energy_ij!(system, i, j, position_i, species_i)
        system.local_energy[j] += energy_ij
        energy += energy_ij
    end
    # Update particle i local energy cache
    system.local_energy[i] = energy
    return energy
end

# With linked list

function destroy_particle!(system::Atoms, i, cell_list::LinkedList)
    energy = zero(typeof(system.density))
    # Get cell of particle i
    position_i, species_i = system.position[i], system.species[i]
    mc = get_cell(position_i, cell_list.cell)
    push!(system.cache, (i, system.local_energy[i]))
    # Scan the neighbourhood of cell mc (including itself)
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(mc2, cell_list.ncells)
        # Scan atoms in cell c2
        j = cell_list.head[c2]
        while (j != -1)
            push!(system.cache, (j, system.local_energy[j]))
            energy_ij = calc_energy_ij!(system, i, j, position_i, species_i)
            system.local_energy[j] -= energy_ij
            energy += energy_ij
            j = cell_list.list[j]
        end
    end
    # Update particle i local energy cache
    system.local_energy[i] = energy
    return energy
end

function create_particle!(system::Atoms, i, cell_list::LinkedList)
    energy = zero(typeof(system.density))
    # Get cell of particle i
    position_i, species_i = system.position[i], system.species[i]
    mc = get_cell(position_i, cell_list.cell)
    # Scan the neighbourhood of cell mc (including itself)
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(mc2, cell_list.ncells)
        # Scan atoms in cell c2
        j = cell_list.head[c2]
        while (j != -1)
            energy_ij = calc_energy_ij!(system, i, j, position_i, species_i)
            system.local_energy[j] += energy_ij
            energy += energy_ij
            j = cell_list.list[j]
        end
    end
    # Update particle i local energy cache
    system.local_energy[i] = energy
    return energy
end

nothing

