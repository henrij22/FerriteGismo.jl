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
    field_index = Ferrite.find_field(only(dh.subdofhandlers), field_name)
    return fixEdge!(dh, ch, side, field_index; components)
end

function fixEdge!(
        dh::IGADofHandler, ch::ConstraintHandler, side::Union{Symbol, Integer},
        field_index::Int; components = nothing
    )
    side_index = isa(side, Symbol) ? _side_idx_for_symbol(side) : side
    field_ip = Ferrite.getfieldinterpolation(dh, (1, field_index))
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


# ==============================================================================
# Inhomogeneous (and time/load dependent) conditions on a parametric edge
# ==============================================================================

"""
    edgeControlPoints(dh::IGADofHandler, side, field_name::Symbol) -> OrderedSet{Int}

Global indices of the control points of `field_name` on one side of the parametric domain.

`side` is one of `:left`, `:right`, `:bottom`, `:top`, or the corresponding G+Smo boundary
index.

These are the "nodes" of an [`IGAGrid`](@ref), i.e. exactly the entities a
[`Dirichlet`](https://ferrite-fem.github.io/Ferrite.jl/stable/reference/boundary_conditions/)
condition on the edge has to prescribe.
"""
function edgeControlPoints(dh::IGADofHandler, side::Union{Symbol, Integer}, field_name::Symbol)
    side_index = isa(side, Symbol) ? _side_idx_for_symbol(side) : side
    field_index = Ferrite.find_field(only(dh.subdofhandlers), field_name)
    field_ip = Ferrite.getfieldinterpolation(dh, (1, field_index))
    base_ip = Ferrite.get_base_interpolation(field_ip)
    return Ferrite.OrderedSet{Int}(
        Int(p) for p in toMatrix(boundary(_maybeDeref(base_ip.basis), side_index))
    )
end

"""
    prescribeEdge!(dh::IGADofHandler, ch::ConstraintHandler, side, field_name::Symbol, f; components = nothing)

Prescribe the (possibly inhomogeneous and time dependent) value `f` for the field
`field_name` on the boundary `side` of the parametric domain — the inhomogeneous
counterpart of [`fixEdge!`](@ref).

`side` is one of `:left`, `:right`, `:bottom`, `:top` (or the corresponding G+Smo boundary
index) and `components` selects the components of a vector-valued field, as for `fixEdge!`.
`f` is either a number or a function `f(x)` / `f(x, t)` of the control-point coordinate `x`
and the time / load factor `t`, returning one value per constrained component. The
condition is registered as a regular Ferrite `Dirichlet` constraint on the control points of
the edge, so `update!(ch, t)` and `apply!`/`apply_zero!` behave as usual, and load-stepping
or path-following code that differentiates the prescribed values with respect to the load
factor keeps working.

!!! note "Control-point values"
    The value is prescribed on the *control points* of the edge, with `f` evaluated at the
    control-point coordinates. Since the spline basis is a partition of unity, an `f` that
    is constant along the edge — the usual displacement-control case — is reproduced
    exactly. A value varying along the edge is only interpolated at the control points, not
    L2-projected onto the edge, so it is met approximately.

# Example
```julia
ch = ConstraintHandler(dh)
# displacement control: u_x = -t on the left edge
prescribeEdge!(dh, ch, :left, :u, (x, t) -> -t; components = 1)
close!(ch)
```
"""
function prescribeEdge!(
        dh::IGADofHandler, ch::ConstraintHandler, side::Union{Symbol, Integer},
        field_name::Symbol, f; components = nothing
    )
    controlPoints = edgeControlPoints(dh, side, field_name)
    ncomp = _nComponents(dh, field_name)
    nprescribed = components === nothing ? ncomp : length(components)
    # A constant is broadcast over the constrained components, so that e.g. `0.0` clamps a
    # vector field just like the function `x -> zero(Vec{ncomp})` would.
    fFunction = if isa(f, Function)
        f
    elseif isa(f, Number) && nprescribed > 1
        (x -> fill(f, nprescribed))
    else
        (x -> f)
    end
    add!(ch, Dirichlet(field_name, controlPoints, fFunction, components))
    return ch
end

function _nComponents(dh::IGADofHandler, field_name::Symbol)
    field_index = Ferrite.find_field(only(dh.subdofhandlers), field_name)
    return Ferrite.n_components(Ferrite.getfieldinterpolation(dh, (1, field_index)))
end

#=
Ferrite's `add!(ch, ::Dirichlet)` derives the constrained dofs from the interpolation of the
cell the entity belongs to, which the spline basis of an IGA field cannot answer. The IGA
"nodes" are the control points though, so a `Dirichlet` on a set of control-point indices
can be resolved directly from the interleaved dof numbering (see the `IGADofHandler`
`close!`). The `Dirichlet` is stored the way Ferrite stores node-set conditions — the
control-point indices in `local_facet_dofs` and the global dofs in `local_facet_dofs_offset`
— so that Ferrite's `update!` path for node sets computes the inhomogeneities.
=#
function Ferrite.add!(ch::ConstraintHandler{<:IGADofHandler}, dbc::Dirichlet)
    @argcheck !Ferrite.isclosed(ch) "The ConstraintHandler is already closed"
    @argcheck eltype(dbc.facets) === Int "Dirichlet conditions on an IGA grid have to be given on a set of control-point indices, e.g. via `prescribeEdge!`"

    dh = ch.dh
    field_index = Ferrite.find_field(only(dh.subdofhandlers), dbc.field_name)
    field_ip = Ferrite.getfieldinterpolation(dh, (1, field_index))
    ncomp = Ferrite.n_components(field_ip)
    offset = dh.field_offsets[field_index]

    isempty(dbc.components) && append!(dbc.components, 1:ncomp)
    @argcheck all(c -> 1 <= c <= ncomp, dbc.components) "Components $(dbc.components) out of range for field :$(dbc.field_name) with $ncomp component(s)"

    constrainedDofs = Int[]
    empty!(dbc.local_facet_dofs)
    for p in dbc.facets
        # Interleaved dof numbering (u1x, u1y, u2x, u2y, ...), as in `fixEdge!`
        for c in dbc.components
            push!(constrainedDofs, offset + ncomp * (p - 1) + c)
        end
        push!(dbc.local_facet_dofs, p)  # control point index, cf. Ferrite's node-set path
    end
    copy!(dbc.local_facet_dofs_offset, constrainedDofs) # global dofs, cf. Ferrite's node-set path

    push!(ch.dbcs, dbc)
    # Unused on the node-set update path, but `ch.bcvalues` is indexed in lockstep with `ch.dbcs`
    push!(ch.bcvalues, eltype(ch.bcvalues)(Array{Float64, 3}(undef, 0, 0, 0), Int[], 0))
    for d in constrainedDofs
        Ferrite.add_prescribed_dof!(ch, d, NaN, nothing)
    end
    return ch
end
