using ConcreteStructs

struct Molecules{D,  VS<:AbstractVector, M<:Model, C<:CellList, T<:AbstractFloat} <: Particles
    position::Vector{SVector{D,T}}
    species::VS
    mol_species::Vector{Tuple{Int, Int, Int, Int}}
    density::T
    temperature::T
    model::M
    d::Int
    N::Int
    N_mol::Int
    box::SVector{D,T}
    local_energy::Vector{T}
    cell_list::C
    cache::Vector{Tuple{Int,T}}
    particle_ids::Base.OneTo{Int}
    mol_ids::Base.OneTo{Int}
    bonds::Vector{Vector{Int}}
end


function System(position, species, mol_species, density::T, temperature::T, model::Model, bonds; list_type=EmptyList) where {T<:AbstractFloat}
    @assert length(position) == length(species)
    particle_ids = eachindex(position)
    mol_ids = eachindex(mol_species)
    N = length(particle_ids)
    N_mol = length(mol_ids)
    d = length(Array(position)[1])
    box = @SVector fill(T((N / density)^(1 / d)), d)
    e_locals = zeros(T, N)
    maxcut = maximum([cutoff(spi, spj, model) for spi in species for spj in species])
    cell_list = list_type(box, maxcut, N)
    cache = Vector{Tuple{Int,T}}()
    system = Molecules(position, species, mol_species, density, temperature, model, d, N, N_mol, box, e_locals, cell_list, cache, particle_ids, mol_ids, bonds)
    build_cell_list!(system, system.cell_list)
    system.local_energy .= [local_energy(system, i, cell_list) for i in particle_ids]
    return system
end

function calc_energy_ij!(system::Molecules, i, j,  position_i, species_i, bonded)
    if i != j
        energy = zero(typeof(system.density))
        position_j, species_j = system.position[j], system.species[j]
        image = nearest_image_distance(position_i, position_j, system.box)
        r2 = dot(image, image)
        if r2 ≤ cutoff2(species_i, species_j, system.model)
            energy += potential(r2, species_i, species_j, system.model)
        end
        if bonded
            energy += bond_potential(r2, species_i, species_j, system.model)
        end
        return energy
    end
    return zero(typeof(system.density))
end

function local_energy(system::Molecules, i, ::EmptyList)
    energy = zero(typeof(system.density))
    position_i, species_i = system.position[i], system.species[i]
    bonds_i = system.bonds[i]
    @inbounds for j in system.particle_ids
        bonded = j ∈ bonds_i
        energy += calc_energy_ij!(system, i, j, position_i, species_i, bonded)
    end
    return energy
end

function local_energy(system::Molecules, i, cell_list::LinkedList)
    energy = zero(typeof(system.density))
    position_i, species_i = system.position[i], system.species[i]
    mc = get_cell(position_i, cell_list.cell)
    # Scan the neighbourhood of cell mc (including itself)
    bonds_i = system.bonds[i]
    for j in bonds_i
        energy += calc_energy_ij!(system, i, j, position_i, species_i, true)
    end
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(mc2, cell_list.ncells)
        # Scan atoms in cell c2
        j = cell_list.head[c2]
        while (j != -1)
            if j ∉ bonds_i
                energy += calc_energy_ij!(system, i, j, position_i, species_i, false)
            end
            j = cell_list.list[j]
        end
    end
    return energy
end

function destroy_particle!(system::Molecules, i, ::EmptyList)
    energy = zero(typeof(system.density))
    position_i, species_i = system.position[i], system.species[i]
    bonds_i = system.bonds[i]
    push!(system.cache, (i, system.local_energy[i]))
    # Loop over particles
    @inbounds for j in system.particle_ids
        push!(system.cache, (j, system.local_energy[j]))
        bonded = j ∈ bonds_i
        energy_ij = calc_energy_ij!(system, i, j, position_i, species_i, bonded)
        system.local_energy[j] -= energy_ij
        energy += energy_ij
    end
    # Update particle i local energy cache
    system.local_energy[i] = energy
    return energy
end

function create_particle!(system::Molecules, i, ::EmptyList)
    energy = zero(typeof(system.density))
    position_i, species_i = system.position[i], system.species[i]
    bonds_i = system.bonds[i]
    # Loop over particles
    @inbounds for j in system.particle_ids
        bonded = j ∈ bonds_i
        energy_ij = calc_energy_ij!(system, i, j, position_i, species_i, bonded)
        system.local_energy[j] += energy_ij
        energy += energy_ij
    end
    # Update particle i local energy cache
    system.local_energy[i] = energy
    return energy
end

# With linked list

function destroy_particle!(system::Molecules, i, cell_list::LinkedList)
    energy = zero(typeof(system.density))
    # Get cell of particle i
    position_i, species_i = system.position[i], system.species[i]
    mc = get_cell(position_i, cell_list.cell)
    push!(system.cache, (i, system.local_energy[i]))
    # Scan the neighbourhood of cell mc (including itself)
    bonds_i = system.bonds[i]
    @inbounds for j in bonds_i
        push!(system.cache, (j, system.local_energy[j]))
        energy_ij = calc_energy_ij!(system, i, j, position_i, species_i, true)
        system.local_energy[j] -= energy_ij
        energy += energy_ij
    end
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(mc2, cell_list.ncells)
        # Scan atoms in cell c2
        j = cell_list.head[c2]
        while (j != -1)
            if j ∉ bonds_i
                push!(system.cache, (j, system.local_energy[j]))
                energy_ij = calc_energy_ij!(system, i, j, position_i, species_i, false)
                system.local_energy[j] -= energy_ij
                energy += energy_ij
            end
            j = cell_list.list[j]
        end
    end
    # Update particle i local energy cache
    system.local_energy[i] = energy
    return energy
end

function create_particle!(system::Molecules, i, cell_list::LinkedList)
    energy = zero(typeof(system.density))
    # Get cell of particle i
    position_i, species_i = system.position[i], system.species[i]
    mc = get_cell(position_i, cell_list.cell)
    # Scan the neighbourhood of cell mc (including itself)
    bonds_i = system.bonds[i]
    @inbounds for j in bonds_i
        energy_ij = calc_energy_ij!(system, i, j, position_i, species_i, true)
        system.local_energy[j] += energy_ij
        energy += energy_ij
    end
    @inbounds for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
        # Calculate the scalar cell index of the neighbour cell (with PBC)
        c2 = cell_index(mc2, cell_list.ncells)
        # Scan atoms in cell c2
        j = cell_list.head[c2]
        while (j != -1)
            if j ∉ bonds_i
                energy_ij = calc_energy_ij!(system, i, j, position_i, species_i, false)
                system.local_energy[j] += energy_ij
                energy += energy_ij
            end
            j = cell_list.list[j]
        end
    end
    # Update particle i local energy cache
    system.local_energy[i] = energy
    return energy
end

nothing

