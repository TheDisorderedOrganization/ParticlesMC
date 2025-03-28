using StaticArrays

abstract type NeighbourList end

function build_cell_list!(system::Particles)
    return build_cell_list!(system, get_neighbour_list(system))
end
struct EmptyList <: NeighbourList end

function EmptyList(box, rcut, N)
    return EmptyList()
end

function build_cell_list!(::Particles,::EmptyList)
    return nothing
end

function update_cell_list!(i, c, c2, ::EmptyList)
    return nothing
end

function old_new_cell(::Particles, i, ::EmptyList)
    return 1, 1
end

struct CellList{T<:AbstractFloat, d} <: NeighbourList
    cs::Vector{Int}
    box::SVector{d,T}
    ncells::NTuple{d,Int}
    cells::Vector{Vector{Int}}  # List of particles in each cell
end

function CellList(box::SVector{d,T}, rcut::T, N::Int) where {d,T<:AbstractFloat}
    
    # Calculate cell dimensions ensuring they're >= rcut
    cell_box = @. box / floor(Int, box / rcut)
    
    # Calculate number of cells in each dimension
    ncells = ntuple(i -> max(1, floor(Int, box[i] / cell_box[i])), d)
    
    # Initialize empty cells
    cells = [Int[] for _ in 1:prod(ncells)]
    
    # Initialize cell indices for particles
    cs = zeros(Int, N)
    
    return CellList(cs, cell_box, ncells, cells)
end

function build_cell_list!(system::Particles, cell_list::CellList)
    # Reset cell list
    #foreach(empty!, cell_list.cells)
    #cell_list.cells = [Int[] for _ in 1:prod(cell_list.ncells)]
    # Populate cell list
    for (i, position_i) in enumerate(system)
        mc = get_cell(position_i, cell_list)
        c = cell_index(cell_list, mc)
        cell_list.cs[i] = c
        push!(cell_list.cells[c], i)  # Directly append particle index
    end
    return nothing
end

function update_cell_list!(i, c, c2, cell_list::CellList)
    # Remove from old cell using in-place filter
    filter!(x -> x != i, cell_list.cells[c])
    # Add to new cell
    push!(cell_list.cells[c2], i)
    # Update particle's cell index
    cell_list.cs[i] = c2
    return nothing
end

function cell_index(cell_list::CellList, mc::NTuple{2,Int})
    i, j = mc          # Cell indices (can be negative)
    ncx, ncy = cell_list.ncells  # Grid dimensions
    # Handle negative indices with proper modulo
    # Julia's mod1 function shifts the result to 1:n range
    i = mod1(i, ncx)
    j = mod1(j, ncy)
    # Convert to linear index (1-based)
    return i + (j-1) * ncx
end

function cell_index(cell_list::CellList, mc::NTuple{3,Int})
    i, j, k = mc          # Cell indices (can be negative)
    ncx, ncy, ncz = cell_list.ncells  # Grid dimensions
    
    # Handle negative indices with mod1
    i = mod1(i, ncx)
    j = mod1(j, ncy)
    k = mod1(k, ncz)
    
    # Convert to linear index (1-based)
    return i + (j-1)*ncx + (k-1)*ncx*ncy
end

function get_cell(xi::SVector{d,T}, box::SVector{d,T}) where {d,T<:AbstractFloat}
    return ntuple(a -> Int(fld(xi[a], box[a])), d)
end

function get_cell(xi::SVector{d,T}, cell_list::NeighbourList) where {d,T<:AbstractFloat}
    return get_cell(xi, cell_list.box)
end


function old_new_cell(system::Particles, i, cell_list::CellList)
    c = system.cell_list.cs[i]
    # New cell index
    mc2 = get_cell(get_position(system, i), cell_list)
    c2 = cell_index(system.cell_list, mc2)
    return c, c2
end

struct LinkedList{T<:AbstractFloat, d} <: NeighbourList
    cs::Vector{Int}
    box::SVector{d,T}
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

function build_cell_list!(system::Particles, cell::LinkedList)
    # Reset the headers
    for c in 1:prod(system.cell_list.ncells)
        system.cell_list.head[c] = -1
    end
    # Scan system to construct headers and linked lists
    for i in eachindex(system)
        # Vector cell indices to which this atom belongs (0 ≤ mc[a] < nc[a] ∀a)
        mc = get_cell(get_position(system, i), system.cell_list.box)
        # Translate the vector cell indices, mc, to a scalar cell index
        c = cell_index(system.cell_list, mc)
        # Link to the previous occupant (or EMPTY (-1) if you're the 1st)
        system.cell_list.list[i] = system.cell_list.head[c]
        # The last one goes to the header
        system.cell_list.head[c] = i
        # Save cell indices
        system.cell_list.cs[i] = c
    end
    return nothing
end


function cell_index(cell_list::LinkedList, mc::NTuple{d,Int}) where {d}
    mc_fold = fold_back(mc, cell_list.ncells)
    c = 1
    stride = 1
    for i in reverse(1:d)
        c += mc_fold[i] * stride
        stride *= cell_list.ncells[i]
    end
    return c
end

function update_cell_list!(i::Int, c::Int, c2::Int, cell_list::LinkedList)
    # Remove particle from old cell (distinguish head or not)
    if cell_list.head[c] == i
        cell_list.head[c] = cell_list.list[i]
    else
        # Find preceding particle in old cell
        j = cell_list.head[c]
        while cell_list.list[j] != i
            j = cell_list.list[j]
        end
        # Update list reference in old cell
        cell_list.list[j] = cell_list.list[i]
    end
    # Insert particle into new cell (exchange with head)
    cell_list.list[i] = cell_list.head[c2]
    cell_list.head[c2] = i
    # Update cell index
    cell_list.cs[i] = c2
    return nothing
end

function old_new_cell(system::Particles, i, cell_list::LinkedList)
    c = system.cell_list.cs[i]
    # New cell index
    mc2 = get_cell(get_position(system, i), cell_list)
    c2 = cell_index(cell_list, mc2)
    return c, c2
end

