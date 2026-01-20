using StaticArrays

abstract type NeighbourList end

"""Build or update the neighbour list for `system` using its current neighbour list type.

This dispatches to the concrete `build_neighbour_list!` method for the active neighbour list implementation.
"""
function build_neighbour_list!(system::Particles)
    return build_neighbour_list!(system, get_neighbour_list(system))
end

"""A no-op neighbour list implementation: always empty.

Useful for testing or systems where no neighbour list is desired.
"""
struct EmptyList <: NeighbourList end

"""Construct an `EmptyList` (ignores box, rcut, N).
"""
function EmptyList(box, rcut, N)
    return EmptyList()
end

"""No-op build for `EmptyList`.
"""
function build_neighbour_list!(::Particles,::EmptyList)
    return nothing
end

"""No-op update for `EmptyList`.
"""
function update_neighbour_list!(i, c, c2, ::EmptyList)
    return nothing
end

"""Return placeholder old and new cell indices for `EmptyList`.

Always returns (1,1) as no cells are tracked.
"""
function old_new_cell(::Particles, i, ::EmptyList)
    return 1, 1
end

"""Calling an EmptyList objects return an object which can be iterated upon.

This iteration will return the indices of the neighbours (which for this list is all the other particles in the system).
"""
function (empty_list::EmptyList)(system::Particles, ::Int)
    return (j for j in 1:length(system))
end


"""Return the scalar cell index of particle `i` stored in `neighbour_list`.
"""
function get_cell_index(i::Int, neighbour_list::NeighbourList)
    return neighbour_list.cs[i]
end
"""Cell-list neighbour list implementation.

Fields:
- `cs`: per-particle scalar cell index
- `box`: cell dimensions
- `ncells`: number of cells per dimension
- `cells`: vectors of particle indices in each cell
- `neighbour_cells`: precomputed neighbouring cell indices for each cell
"""
struct CellList{T<:AbstractFloat, d} <: NeighbourList
    cs::Vector{Int}
    box::SVector{d,T}
    ncells::NTuple{d,Int}
    cells::Vector{Vector{Int}}  # List of particles in each cell
    neighbour_cells::Vector{Vector{Int}}  # List of neighbouring cells
end

"""Return the scalar cell index corresponding to `mc` using `neighbour_list` metadata.
"""
function cell_index(neighbour_list::NeighbourList, mc::NTuple{d,Int}) where {d}
    return cell_index(neighbour_list.ncells, mc)
end

"""Compute a scalar cell index from the multi-dimensional cell coordinate `mc`.

Uses row-major ordering with periodic boundary treatment via `fold_back`.
"""
function cell_index(ncells::NTuple{d, Int}, mc::NTuple{d,Int}) where {d}
    mc_fold = fold_back(mc, ncells)
    c = 1
    stride = 1
    for i in reverse(1:d)
        c += mc_fold[i] * stride
        stride *= ncells[i]
    end
    return c
end

"""Build for each scalar cell index the list of neighbouring cells (including itself).

Returns a vector where element `c` contains a vector of scalar indices of neighbouring cells.
"""
function build_neighbour_cells(ncells::NTuple{d,Int}) where {d}
    # Create a list of neighbour cells for each cell
    neighbourcells = Vector{Vector{Int}}(undef, prod(ncells))
    ranges = map(Base.OneTo, ncells)  # create ranges 1:n1, 1:n2, ...
    for mc in Iterators.product(ranges...)
        c = cell_index(ncells, mc)
        neighbours = []
        for mc2 in Iterators.product(map(x -> x-1:x+1, mc)...)
            # Calculate the scalar cell index of the neighbour cell (with PBC)
            c2 = cell_index(ncells, mc2)
            push!(neighbours, c2)
        end
        # The unique is required if there are fewer than 3 cells in one direction
        # Otherwise, the same index can appear multiple times
        neighbourcells[c] = unique!(neighbours)
    end
    return neighbourcells
end

"""Construct a `CellList` neighbour list given a box, cutoff `rcut`, and particle count `N`.

Cells are chosen so that their dimensions are at least `rcut`, and neighbour cells are precomputed.
"""
function CellList(box::SVector{d,T}, rcut::T, N::Int) where {d,T<:AbstractFloat}

    # Calculate cell dimensions ensuring they're >= rcut
    cell_box = @. box / floor(Int, box / rcut)

    # Calculate number of cells in each dimension
    ncells = ntuple(i -> max(1, floor(Int, box[i] / cell_box[i])), d)

    # Initialize empty cells
    cells = [Int[] for _ in 1:prod(ncells)]

    # Initialize cell indices for particles
    cs = zeros(Int, N)
    neighbour_cells = build_neighbour_cells(ncells)
    return CellList(cs, cell_box, ncells, cells, neighbour_cells)
end

"""Populate `neighbour_list` (a `CellList`) with particle indices from `system`.

Assigns each particle to its corresponding cell and appends it to that cell's list.
"""
function build_neighbour_list!(system::Particles, neighbour_list::CellList)
    # Reset cell list
    #foreach(empty!, neighbour_list.cells)
    #neighbour_list.cells = [Int[] for _ in 1:prod(neighbour_list.ncells)]
    # Populate cell list
    for (i, position_i) in enumerate(system)
        c = get_cell_index(position_i, neighbour_list)
        neighbour_list.cs[i] = c
        push!(neighbour_list.cells[c], i)  # Directly append particle index
    end
    return nothing
end

"""Update particle `i` moving from cell `c` to cell `c2` in a `CellList`.

Performs an in-place removal from the old cell and appends to the new one.
"""
function update_neighbour_list!(i, c, c2, neighbour_list::CellList)
    # Remove from old cell using in-place filter
    filter!(x -> x != i, neighbour_list.cells[c])
    # Add to new cell
    push!(neighbour_list.cells[c2], i)
    # Update particle's cell index
    neighbour_list.cs[i] = c2
    return nothing
end




"""Return the multi-dimensional cell coordinate for position `xi` in `box`.

Result is an `NTuple{d,Int}` giving the cell index in each dimension.
"""
function get_cell(xi::SVector{d,T}, box::SVector{d,T}) where {d,T<:AbstractFloat}
    return ntuple(a -> Int(fld(xi[a], box[a])), d)
end

"""Return the multi-dimensional cell coordinate for position `xi` using `neighbour_list` metadata.
"""
function get_cell(xi::SVector{d,T}, neighbour_list::NeighbourList) where {d,T<:AbstractFloat}
    return get_cell(xi, neighbour_list.box)
end

"""Return the scalar cell index for position `xi` using `neighbour_list`.
"""
function get_cell_index(xi::SVector{d,T}, neighbour_list::NeighbourList) where {d,T<:AbstractFloat}
    mc = get_cell(xi, neighbour_list)
    return cell_index(neighbour_list, mc)
end
"""Return the old and new cell indices for particle `i` in `neighbour_list` (CellList).

Used to detect whether a particle has moved between cells.
"""
function old_new_cell(system::Particles, i, neighbour_list::CellList)
    c = neighbour_list.cs[i]
    # New cell index
    mc2 = get_cell(get_position(system, i), neighbour_list)
    c2 = cell_index(neighbour_list, mc2)
    return c, c2
end

"""Calling an EmptyList objects return an object which can be iterated upon.

This iteration will return the indices of the neighbours (which for this list is all the other particles in the system).
"""
function (cell_list::CellList)(system::Particles, i::Int)
    position_i = get_position(system, i)
    c = get_cell_index(position_i, cell_list)
    neighbour_cells = cell_list.neighbour_cells[c]
    # Scan the neighbourhood of cell mc (including itself)
    # and from there scan atoms in cell c2
    return (j for c2 in neighbour_cells for j in @inbounds cell_list.cells[c2])
end


"""Linked-list neighbour list implementation.

Uses arrays `head` and `list` to store per-cell linked lists of particle indices.
"""
struct LinkedList{T<:AbstractFloat, d} <: NeighbourList
    cs::Vector{Int}
    box::SVector{d,T}
    ncells::NTuple{d,Int}
    head::Vector{Int}
    list::Vector{Int}
    rcut::T
    neighbour_cells::Vector{Vector{Int}}
end

"""Construct a `LinkedList` neighbour list given box, cutoff `rcut`, and number of particles `N`.
"""
function LinkedList(box, rcut, N)
    cell = box ./ fld.(box, rcut)
    ncells = ntuple(a -> Int(box[a] / cell[a]), length(box))
    head = -ones(Int, prod(ncells))
    list = zeros(Int, N)
    cs = zeros(Int, N)
    neighbour_cells = build_neighbour_cells(ncells)
    return LinkedList(cs, cell, ncells, head, list, rcut, neighbour_cells)
end


"""Populate `neighbour_list` (a `LinkedList`) by constructing head/list linked lists per cell.

Resets headers and inserts particles into per-cell linked lists efficiently.
"""
function build_neighbour_list!(system::Particles, neighbour_list::LinkedList)
    # Reset the headers
    for c in 1:prod(neighbour_list.ncells)
        neighbour_list.head[c] = -1
    end
    # Scan system to construct headers and linked lists
    for i in eachindex(system)
        # Vector cell indices to which this atom belongs (0 ≤ mc[a] < nc[a] ∀a)
        mc = get_cell(get_position(system, i), neighbour_list.box)
        # Translate the vector cell indices, mc, to a scalar cell index
        c = cell_index(neighbour_list, mc)
        # Link to the previous occupant (or EMPTY (-1) if you're the 1st)
        neighbour_list.list[i] = neighbour_list.head[c]
        # The last one goes to the header
        neighbour_list.head[c] = i
        # Save cell indices
        neighbour_list.cs[i] = c
    end
    return nothing
end

"""Update particle `i` moving from cell `c` to `c2` in a `LinkedList`.

Performs linked-list updates handling header and interior-node cases.
"""
function update_neighbour_list!(i::Int, c::Int, c2::Int, neighbour_list::LinkedList)
    # Remove particle from old cell (distinguish head or not)
    if neighbour_list.head[c] == i
        neighbour_list.head[c] = neighbour_list.list[i]
    else
        # Find preceding particle in old cell
        j = neighbour_list.head[c]
        while neighbour_list.list[j] != i
            j = neighbour_list.list[j]
        end
        # Update list reference in old cell
        neighbour_list.list[j] = neighbour_list.list[i]
    end
    # Insert particle into new cell (exchange with head)
    neighbour_list.list[i] = neighbour_list.head[c2]
    neighbour_list.head[c2] = i
    # Update cell index
    neighbour_list.cs[i] = c2
    return nothing
end

"""Return old and new cell indices for particle `i` using a `LinkedList` neighbour list.
"""
function old_new_cell(system::Particles, i, neighbour_list::LinkedList)
    c = neighbour_list.cs[i]
    # New cell index
    mc2 = get_cell(get_position(system, i), neighbour_list)
    c2 = cell_index(neighbour_list, mc2)
    return c, c2
end

""" This struct is used to iterate over neighbours of a Linked list
"""
struct LinkedIterator
    neighbour_cells::Vector{Int}
    head::Vector{Int}
    list::Vector{Int}
end

# To iterate over the neighbours of a linked list, one could write the following loops
#@inbounds for c in neighbour_list.neighbour_cells
#    j = neighbour_list.head[c]
#    while (j != -1)
#        do stuff
#        j = neighbour_list.list[j]
#    end
#end
# This is however impossible to rewrite as a simple generator
# So we implement the following function, which uses a state to carry over the needed information
function Base.iterate(neighbour_list::LinkedIterator, state=-1)
    # First time in
    if state == -1
        next = iterate(neighbour_list.neighbour_cells)
        if next == nothing
            return nothing
        end
        c, c_state = next
        j = neighbour_list.head[c]
    else
        c_state, j = state
        j = neighbour_list.list[j]
        if j == -1
            next = iterate(neighbour_list.neighbour_cells, c_state)
            if next == nothing
                return nothing
            end
            c, c_state = next
            j = neighbour_list.head[c]
        end
    end

    if j == -1
        return nothing
    end
    state = (c_state, j)
    return j, state
end


"""Iterate over the particles from adjacent cells.
"""
function get_neighbour_indices(system::Particles, neighbour_list::LinkedList, i::Int)
    position_i = get_position(system, i)
    c = get_cell_index(position_i, neighbour_list)
    neighbour_cells = neighbour_list.neighbour_cells[c]

    iterator = LinkedIterator(neighbour_cells, neighbour_list.head, neighbour_list.list)
    return (j for j in iterator)

end
