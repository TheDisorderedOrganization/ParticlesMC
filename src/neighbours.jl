using StaticArrays

abstract type NeighbourList end

function build_neighbour_list!(system::Particles)
    return build_neighbour_list!(system, get_neighbour_list(system))
end

struct EmptyList <: NeighbourList end

function EmptyList(box, rcut, N)
    return EmptyList()
end

function build_neighbour_list!(::Particles,::EmptyList)
    return nothing
end

function update_neighbour_list!(i, c, c2, ::EmptyList)
    return nothing
end

function old_new_cell(::Particles, i, ::EmptyList)
    return 1, 1
end

function get_cell_index(i::Int, neighbour_list::NeighbourList)
    return neighbour_list.cs[i]
end
struct CellList{T<:AbstractFloat, d} <: NeighbourList
    cs::Vector{Int}
    box::SVector{d,T}
    ncells::NTuple{d,Int}
    cells::Vector{Vector{Int}}  # List of particles in each cell
    neighbour_cells::Vector{Vector{Int}}  # List of neighbouring cells
end

function cell_index(neighbour_list::NeighbourList, mc::NTuple{d,Int}) where {d}
    return cell_index(neighbour_list.ncells, mc)
end

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
        neighbourcells[c] = neighbours
    end
    return neighbourcells
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
    neighbour_cells = build_neighbour_cells(ncells)
    return CellList(cs, cell_box, ncells, cells, neighbour_cells)
end

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

function update_neighbour_list!(i, c, c2, neighbour_list::CellList)
    # Remove from old cell using in-place filter
    filter!(x -> x != i, neighbour_list.cells[c])
    # Add to new cell
    push!(neighbour_list.cells[c2], i)
    # Update particle's cell index
    neighbour_list.cs[i] = c2
    return nothing
end




function get_cell(xi::SVector{d,T}, box::SVector{d,T}) where {d,T<:AbstractFloat}
    return ntuple(a -> Int(fld(xi[a], box[a])), d)
end

function get_cell(xi::SVector{d,T}, neighbour_list::NeighbourList) where {d,T<:AbstractFloat}
    return get_cell(xi, neighbour_list.box)
end

function get_cell_index(xi::SVector{d,T}, neighbour_list::NeighbourList) where {d,T<:AbstractFloat}
    mc = get_cell(xi, neighbour_list)
    return cell_index(neighbour_list, mc)
end
function old_new_cell(system::Particles, i, neighbour_list::CellList)
    c = neighbour_list.cs[i]
    # New cell index
    mc2 = get_cell(get_position(system, i), neighbour_list)
    c2 = cell_index(neighbour_list, mc2)
    return c, c2
end

struct LinkedList{T<:AbstractFloat, d} <: NeighbourList
    cs::Vector{Int}
    box::SVector{d,T}
    ncells::NTuple{d,Int}
    head::Vector{Int}
    list::Vector{Int}
    rcut::T
    neighbour_cells::Vector{Vector{Int}}
end

function LinkedList(box, rcut, N)
    cell = box ./ fld.(box, rcut)
    ncells = ntuple(a -> Int(box[a] / cell[a]), length(box))
    head = -ones(Int, prod(ncells))
    list = zeros(Int, N)
    cs = zeros(Int, N)
    neighbour_cells = build_neighbour_cells(ncells)
    return LinkedList(cs, cell, ncells, head, list, rcut, neighbour_cells)
end


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

function old_new_cell(system::Particles, i, neighbour_list::LinkedList)
    c = neighbour_list.cs[i]
    # New cell index
    mc2 = get_cell(get_position(system, i), neighbour_list)
    c2 = cell_index(neighbour_list, mc2)
    return c, c2
end

