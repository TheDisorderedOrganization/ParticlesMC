"""
External fields acting on molecular observables.

Field energies are kept separate from pair-interaction energies because a
molecular field contribution is counted once per molecule, whereas the sum of
per-particle pair energies is divided by two.
"""

using LinearAlgebra
using StaticArrays

"""Abstract type for an external-field contribution to the Hamiltonian."""
abstract type ExternalField end

"""
Marker used when a molecular system has no external field.

Using a concrete no-op field instead of `nothing` keeps field-energy calls
branch-free and preserves the behavior of systems constructed without a field.
"""
struct NoExternalField <: ExternalField end

"""
    AligningField(h, orientation_definition)

An aligning field whose contribution for molecule `m` is

```math
U_m = -\\mathbf{h} \\cdot \\mathbf{n}_m,
```

where `n_m` is the unit vector computed using `orientation_definition`. This linear coupling is polar: reversing `n_m` changes the
sign of the field energy.
"""
struct AligningField{H<:SVector,O<:OrientationDefinition} <: ExternalField
    h::H
    orientation_definition::O
end

function AligningField(h::AbstractVector{<:Real},
                       orientation_definition::OrientationDefinition)
    isempty(h) && throw(ArgumentError("the external-field vector cannot be empty"))
    all(isfinite, h) || throw(ArgumentError("external-field components must be finite"))
    h_float = float.(h)
    return AligningField(SVector{length(h_float)}(h_float), orientation_definition)
end

"""
    field_energy(field, system, molecule_id)

Return the external-field contribution of one molecule. Move implementations
use this local quantity to compute a field-energy difference without summing
over unaffected molecules.
"""
field_energy(::NoExternalField, system, molecule_id::Int) = zero(typeof(system.temperature))

function field_energy(field::AligningField, system, molecule_id::Int)
    n = orientation(field.orientation_definition, system, molecule_id)
    return -dot(field.h, n)
end

"""
    total_field_energy(field, system)

Return the full one-body field contribution

```math
U_{\\mathrm{field}} = \\sum_{m=1}^{N_{\\mathrm{mol}}} U_m.
```

Each molecular term is counted exactly once. In particular, this result must
not be included in the per-particle pair-energy sum that is later halved.
"""
function total_field_energy(field::ExternalField, system)
    return mapreduce(
        molecule_id -> field_energy(field, system, molecule_id),
        +,
        Base.OneTo(system.Nmol);
        init=zero(typeof(system.temperature)),
    )
end

"""Return the field-vector dimension, or `nothing` for a dimensionless no-op field."""
field_dimension(::NoExternalField) = nothing
field_dimension(field::AligningField) = length(field.h)

"""Convert an absent or already constructed field specification to a field object."""
external_field_from_config(::Nothing) = NoExternalField()
external_field_from_config(field::ExternalField) = field

"""
Build an external field from a TOML-compatible dictionary. This is the boundary
between user-facing configuration values and the typed field/orientation
objects used by the simulation.

Supported orientations are `"plane_normal"` with `atoms = [a,b,c]` and
`"center_to_atom"` with `atom = a`.
"""
function external_field_from_config(config::AbstractDict)
    field_type = get(config, "type", "align")
    field_type == "align" || throw(ArgumentError("unsupported external field type: $field_type"))
    haskey(config, "h") || throw(ArgumentError("an aligning field requires an `h` vector"))

    orientation_name = get(config, "orientation", "plane_normal")
    definition = if orientation_name == "plane_normal"
        atoms = get(config, "atoms", [1, 2, 3])
        length(atoms) == 3 || throw(ArgumentError("plane-normal orientation requires three atoms"))
        PlaneNormalOrientation(atoms...)
    elseif orientation_name == "center_to_atom"
        CenterToAtomOrientation(get(config, "atom", 1))
    else
        throw(ArgumentError("unsupported orientation definition: $orientation_name"))
    end

    return AligningField(config["h"], definition)
end
