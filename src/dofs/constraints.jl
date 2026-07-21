function _side_idx_for_symbol(side::Symbol)::Int
    if side == :left
        return 1
    elseif side == :right
        return 2
    elseif side == :bottom
        return 3
    elseif side == :top
        return 4
    else
        throw(ArgumentError("Unknown side $side"))
    end
end


"""
    fixEdge!(dh::IGADofHandler, ch::ConstraintHandler, side, field_name::Symbol; components = nothing)
    fixEdge!(dh::IGADofHandler, ch::ConstraintHandler, side, field_index::Int; components = nothing)

Prescribe a homogeneous (zero) Dirichlet condition on the degrees of freedom of the field
`field_name` (or `field_index`) that are associated with the boundary `side` of the
parametric domain.

`side` is either one of the symbols `:left`, `:right`, `:bottom`, `:top` or the
corresponding G+Smo boundary index. For vector-valued fields, `components` selects which
components to constrain (e.g. `components = 1` to only fix the ``x``-component); by default
all components are constrained. The prescribed dofs are added directly to the constraint
handler `ch`; call `close!(ch)` afterwards as usual.

!!! note
    This is a lightweight helper that works around the fact that the standard Ferrite
    `ConstraintHandler` / `Dirichlet` machinery is not yet fully supported for IGA grids.
    Only homogeneous conditions along whole parametric edges are supported.

# Example
```julia
ch = ConstraintHandler(dh)
fixEdge!(dh, ch, :left, :u)                 # scalar field, whole left edge
fixEdge!(dh, ch, :left, :u; components = 1) # only the x-component of a vector field
close!(ch)
```
"""
function fixEdge!(
        dh::IGADofHandler, ch::ConstraintHandler, side::Union{Symbol, Integer},
        field_name::Symbol; components = nothing
    )
    field_index = Ferrite.find_field(only(dh.dh.subdofhandlers), field_name)
    return fixEdge!(dh, ch, side, field_index; components)
end

function fixEdge!(
        dh::IGADofHandler, ch::ConstraintHandler, side::Union{Symbol, Integer},
        field_index::Int; components = nothing
    )
    side_index = isa(side, Symbol) ? _side_idx_for_symbol(side) : side
    field_ip = Ferrite.getfieldinterpolation(dh.dh, (1, field_index))
    ncomp = Ferrite.n_components(field_ip)
    base_ip = Ferrite.get_base_interpolation(field_ip)
    comps = components === nothing ? (1:ncomp) : components
    offset = dh.field_offsets[field_index]

    boundary_cps = toMatrix(boundary(_maybeDeref(base_ip.basis), side_index))
    for p in boundary_cps
        # Global control-point index of `p` within this field, matching the interleaved
        # dof numbering (u1x, u1y, u2x, u2y, ...) used when closing the IGADofHandler.
        for c in comps
            dof = offset + ncomp * (p - 1) + c
            Ferrite.add_prescribed_dof!(ch, dof, 0.0)
        end
    end
    return nothing
end
