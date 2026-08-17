"""
Definitions and geometry helpers for molecular orientation vectors.

Orientation definitions are deliberately independent of field energies and Monte
Carlo moves so that the geometrical definition can be changed without modifying
either of those components.
"""

using LinearAlgebra
using StaticArrays

"""Abstract type for a molecular orientation convention."""
abstract type OrientationDefinition end

"""
    CenterToAtomOrientation(atom=1)

Define the molecular orientation as the unit vector from the molecule's
centre of mass to local atom `atom`. The atom index is local to the molecule:
`atom=1` selects the first atom of that molecule, not particle 1 of the system.
"""
struct CenterToAtomOrientation <: OrientationDefinition
    atom::Int

    function CenterToAtomOrientation(atom::Int=1)
        atom > 0 || throw(ArgumentError("the local atom index must be positive"))
        new(atom)
    end
end

"""
    PlaneNormalOrientation(atom1=1, atom2=2, atom3=3)

Define the orientation of a triangular molecule by the right-handed unit normal
to the ordered atoms `(atom1, atom2, atom3)`. Exchanging any two atom indices
reverses the orientation.
"""
struct PlaneNormalOrientation <: OrientationDefinition
    atom1::Int
    atom2::Int
    atom3::Int

    function PlaneNormalOrientation(atom1::Int=1, atom2::Int=2, atom3::Int=3)
        all(i -> i > 0, (atom1, atom2, atom3)) ||
            throw(ArgumentError("local atom indices must be positive"))
        length(unique((atom1, atom2, atom3))) == 3 ||
            throw(ArgumentError("the plane must be defined by three distinct atoms"))
        new(atom1, atom2, atom3)
    end
end

function _normalise_orientation(v)
    magnitude = norm(v)
    magnitude > sqrt(eps(eltype(v))) ||
        throw(ArgumentError("cannot define an orientation from degenerate molecular coordinates"))
    return v / magnitude
end

"""
    unwrap_molecule(positions, box)

Return mutually consistent periodic images of a molecule's positions. The first
atom is the reference and every other atom is placed in its nearest image.
This assumes each molecule is smaller than half the box length in each relevant
direction.
"""
function unwrap_molecule(positions::AbstractVector{<:SVector{D,T}},
                         box::SVector{D,T}) where {D,T<:AbstractFloat}
    isempty(positions) && throw(ArgumentError("a molecule must contain at least one atom"))
    reference = first(positions)
    return [reference + vector(position, reference, box) for position in positions]
end

"""
    molecule_center_of_mass(positions, masses, box)

Compute a molecule's mass-weighted centre, accounting for periodic boundaries.
The returned centre of mass is in the same unwrapped image as the first atom.
"""
function molecule_center_of_mass(positions::AbstractVector{<:SVector{D,T}},
                                 masses::AbstractVector{<:Real},
                                 box::SVector{D,T}) where {D,T<:AbstractFloat}
    length(masses) == length(positions) ||
        throw(ArgumentError("one mass per molecular position is required"))
    all(m -> isfinite(m) && m > zero(m), masses) ||
        throw(ArgumentError("particle masses must be finite and strictly positive"))
    unwrapped = unwrap_molecule(positions, box)
    total_mass = sum(masses)
    return sum(mass * position for (mass, position) in zip(masses, unwrapped)) / total_mass
end

"""
    molecule_center_of_mass(system, molecule_id)

Compute the periodic center of mass of molecule `molecule_id`. This is the
system-level convenience overload of `molecule_center_of_mass(positions,
masses, box)`.
"""
function molecule_center_of_mass(system, molecule_id::Int)
    1 <= molecule_id <= system.Nmol ||
        throw(BoundsError(Base.OneTo(system.Nmol), molecule_id))
    first_atom, last_atom = get_start_end_mol(system, molecule_id)
    return molecule_center_of_mass(
        @view(system.position[first_atom:last_atom]),
        @view(system.mass[first_atom:last_atom]),
        system.box,
    )
end

"""
Compute a center-to-atom orientation from positions and masses belonging to one
molecule. Atom indices in `definition` are local to that molecule.
"""
function orientation(definition::CenterToAtomOrientation,
                     positions::AbstractVector{<:SVector{D,T}},
                     masses::AbstractVector{<:Real},
                     box::SVector{D,T}) where {D,T<:AbstractFloat}
    definition.atom <= length(positions) ||
        throw(BoundsError(positions, definition.atom))
    unwrapped = unwrap_molecule(positions, box)
    center = molecule_center_of_mass(positions, masses, box)
    return _normalise_orientation(unwrapped[definition.atom] - center)
end

"""
Compute a plane-normal orientation from positions belonging to one molecule.
Atom indices in `definition` are local to that molecule.
"""
function orientation(definition::PlaneNormalOrientation,
                     positions::AbstractVector{<:SVector{3,T}},
                     box::SVector{3,T}) where {T<:AbstractFloat}
    indices = (definition.atom1, definition.atom2, definition.atom3)
    maximum(indices) <= length(positions) ||
        throw(BoundsError(positions, maximum(indices)))
    unwrapped = unwrap_molecule(positions, box)
    r1, r2, r3 = (unwrapped[i] for i in indices)
    return _normalise_orientation(cross(r2 - r1, r3 - r1))
end

"""
    orientation(definition, system, molecule_id)

Compute the orientation of molecule `molecule_id` in a `Molecules` system.
Atom indices stored in an orientation definition are local to that molecule.
"""
function orientation(definition::OrientationDefinition,
                     system, molecule_id::Int)
    1 <= molecule_id <= system.Nmol ||
        throw(BoundsError(Base.OneTo(system.Nmol), molecule_id))
    first_atom, last_atom = get_start_end_mol(system, molecule_id)
    return orientation(definition, @view(system.position[first_atom:last_atom]), system.box)
end

function orientation(definition::CenterToAtomOrientation,
                     system, molecule_id::Int)
    1 <= molecule_id <= system.Nmol ||
        throw(BoundsError(Base.OneTo(system.Nmol), molecule_id))
    first_atom, last_atom = get_start_end_mol(system, molecule_id)
    return orientation(
        definition,
        @view(system.position[first_atom:last_atom]),
        @view(system.mass[first_atom:last_atom]),
        system.box,
    )
end
