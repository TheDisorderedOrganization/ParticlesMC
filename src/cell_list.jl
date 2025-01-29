using StaticArrays

abstract type CellList end

struct EmptyList <: CellList end

function EmptyList(box, rcut, N)
    return EmptyList()
end

struct LinkedList{T<:AbstractFloat,d} <: CellList
    cs::Vector{Int}
    cell::SVector{d,T}
    ncells::NTuple{d,Int}
    head::Vector{Int}
    list::Vector{Int}
    rcut::T
end

function LinkedList(box, rcut, N)
    cell = box ./ fld.(box, rcut)
    ncells = ntuple(a -> Int(box[a] / cell[a]), length(box))
    head = -ones(Int, prod(ncells))
    list = zeros(Int, N)
    cs = zeros(Int, N)
    return LinkedList(cs, cell, ncells, head, list, rcut)
end

fold_back(x, box) = x .- fld.(x, box) .* box

function get_cell(xi::SVector{d,T}, cell::SVector{d,T}) where {d,T<:AbstractFloat}
    return ntuple(a -> Int(fld(xi[a], cell[a])), d)
end

function cell_index(mc::NTuple{d,Int}, ncells::NTuple{d,Int}) where {d}
    mc_fold = fold_back(mc, ncells)
    c = 1
    stride = 1
    for i in reverse(1:d)
        c += mc_fold[i] * stride
        stride *= ncells[i]
    end
    return c
end

function build_cell_list!(::Particles, ::EmptyList)
    return nothing
end

function build_cell_list!(system::Particles, ::LinkedList)
    # Reset the headers
    for c in 1:prod(system.cell_list.ncells)
        system.cell_list.head[c] = -1
    end
    # Scan system to construct headers and linked lists
    for i in system.particle_ids
        # Vector cell indices to which this atom belongs (0 ≤ mc[a] < nc[a] ∀a)
        mc = get_cell(system.position[i], system.cell_list.cell)
        # Translate the vector cell indices, mc, to a scalar cell index
        c = cell_index(mc, system.cell_list.ncells)
        # Link to the previous occupant (or EMPTY (-1) if you're the 1st)
        system.cell_list.list[i] = system.cell_list.head[c]
        # The last one goes to the header
        system.cell_list.head[c] = i
        # Save cell indices
        system.cell_list.cs[i] = c
    end
    return nothing
end

function update_cell_list!(::Particles, i, ::EmptyList)
    return nothing
end

function update_cell_list!(system::Particles, i, ::LinkedList)
    # Old cell index
    c = system.cell_list.cs[i]
    # New cell index
    mc2 = get_cell(system.position[i], system.cell_list.cell)
    c2 = cell_index(mc2, system.cell_list.ncells)
    # Check if particle changed cell, otherwise do nothing
    if c != c2
        # Remove particle from old cell (distinguish head or not)
        if system.cell_list.head[c] == i
            system.cell_list.head[c] = system.cell_list.list[i]
        else
            # Find preceding particle in old cell
            j = system.cell_list.head[c]
            while system.cell_list.list[j] != i
                j = system.cell_list.list[j]
            end
            # Update list reference in old cell
            system.cell_list.list[j] = system.cell_list.list[i]
        end
        # Insert particle into new cell (exchange with head)
        system.cell_list.list[i] = system.cell_list.head[c2]
        system.cell_list.head[c2] = i
        # Update cell index
        system.cell_list.cs[i] = c2
    end
    return nothing
end

nothing