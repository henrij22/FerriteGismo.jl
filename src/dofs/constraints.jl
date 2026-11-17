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
    return fixEdge!(dh, ch, side, _globalFieldIndex(dh, field_name); components)
end

function fixEdge!(
        dh::IGADofHandler, ch::ConstraintHandler, side::Union{Symbol, Integer},
        field_index::Int; components = nothing
    )
    side_index = isa(side, Symbol) ? _side_idx_for_symbol(side) : side
    field_ip = _fieldInterpolation(dh, dh.field_names[field_index])
    ncomp = Ferrite.n_components(field_ip)
    base_ip = Ferrite.get_base_interpolation(field_ip)
    comps = components === nothing ? (1:ncomp) : components
    offset = dh.field_offsets[field_index]

    boundary_cps = toMatrix(boundary(_maybeDeref(base_ip.basis), side_index))
    for p in boundary_cps
        # Interleaved dof numbering u1x, u1y, u2x, u2y, ...
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
    base_ip = Ferrite.get_base_interpolation(_fieldInterpolation(dh, field_name))
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
    return Ferrite.n_components(_fieldInterpolation(dh, field_name))
end

#=
Ferrite derives the constrained dofs from the cell interpolation, which a spline basis cannot
answer. The dofs of a control-point set are known directly though, and the `Dirichlet` is
stored like a Ferrite node-set condition (control points in `local_facet_dofs`, global dofs in
`local_facet_dofs_offset`) so that Ferrite's `update!` computes the inhomogeneities.
=#
function Ferrite.add!(ch::ConstraintHandler{<:IGADofHandler, Tv, Ti}, dbc::Dirichlet) where {Tv, Ti}
    @argcheck !Ferrite.isclosed(ch) "The ConstraintHandler is already closed"
    @argcheck eltype(dbc.facets) === Int "Dirichlet conditions on an IGA grid have to be given on a set of control-point indices, e.g. via `prescribeEdge!`"

    dh = ch.dh
    field_ip = _fieldInterpolation(dh, dbc.field_name)
    ncomp = Ferrite.n_components(field_ip)
    offset = fieldOffset(dh, dbc.field_name)

    isempty(dbc.components) && append!(dbc.components, 1:ncomp)
    @argcheck all(c -> 1 <= c <= ncomp, dbc.components) "Components $(dbc.components) out of range for field :$(dbc.field_name) with $ncomp component(s)"

    controlPoints = dbc.facets::Ferrite.OrderedSet{Int}
    constrainedDofs = Int[]
    empty!(dbc.local_facet_dofs)
    for p in controlPoints
        # Interleaved dof numbering (u1x, u1y, u2x, u2y, ...), as in `fixEdge!`
        for c in dbc.components
            push!(constrainedDofs, offset + ncomp * (p - 1) + c)
        end
        push!(dbc.local_facet_dofs, p)
    end
    copy!(dbc.local_facet_dofs_offset, constrainedDofs)

    push!(ch.dbcs, dbc)
    # Unused on the node-set update path, but `ch.bcvalues` is indexed in lockstep with `ch.dbcs`
    push!(ch.bcvalues, Ferrite.BCValues{Tv, Ti}(Array{Tv, 3}(undef, 0, 0, 0), Ti[], zero(Ti)))
    for d in constrainedDofs
        Ferrite.add_prescribed_dof!(ch, d, NaN, nothing)
    end
    return ch
end

# ==============================================================================
# Tangent fixity (clamped edges)
# ==============================================================================

"""
    fixEdgeTangent!(dh::IGADofHandler, ch::ConstraintHandler, side, field_name::Symbol; components = nothing)
    fixEdgeTangent!(dh::IGADofHandler, ch::ConstraintHandler, side, fields::Tuple)

Fix the slope across the boundary `side`, e.g. to clamp plates, beams or Kirchhoff-Love shells.

The `Symbol` form fixes the control points next to `side` to zero, optionally only the given
`components`. Together with `fixEdge!` on
the same side this clamps the edge.

The `Tuple` form is for shell patches (`IGAGrid{3, 2}`), with one 3D vector field `(:u,)` or
three scalar fields `(:ux, :uy, :uz)`. Each control point `d₁` next to `side` gets
`n ⋅ (d₁ - d₀) = 0`, with `d₀` its neighbour on the side and `n` the surface normal there, so
only the normal displacement is tied to the edge. The edge itself may be fixed, prescribed or
free. These are affine constraints, so allocate the matrix with `allocate_matrix(dh, ch)`.

Where two clamped sides meet, both edges have to be held and their prescribed motions must agree
along the normal, otherwise only one of the two conditions is enforced.

Hierarchical (THB) patches are not supported.

# Example
```julia
ch = ConstraintHandler(dh)
fixEdge!(dh, ch, :left, :u)
fixEdgeTangent!(dh, ch, :left, (:u,))
close!(ch)
K = allocate_matrix(dh, ch)
```
"""
function fixEdgeTangent!(
        dh::IGADofHandler, ch::ConstraintHandler, side::Union{Symbol, Integer},
        field_name::Symbol; components = nothing
    )
    ncomp = _nComponents(dh, field_name)
    comps = components === nothing ? (1:ncomp) : components
    @argcheck all(c -> 1 <= c <= ncomp, comps) "Components $comps out of range for field :$field_name with $ncomp component(s)"
    offset = fieldOffset(dh, field_name)
    _, inward = _edgeAndInwardControlPoints(dh, side, field_name)
    for p in inward, c in comps
        Ferrite.add_prescribed_dof!(ch, offset + ncomp * (p - 1) + c, 0.0)
    end
    return nothing
end

function fixEdgeTangent!(
        dh::IGADofHandler, ch::ConstraintHandler, side::Union{Symbol, Integer},
        fields::Tuple{Vararg{Symbol}}
    )
    @argcheck Ferrite.get_grid(dh) isa IGAGrid{3, 2} "the tuple form of fixEdgeTangent! is for shell patches, IGAGrid{3, 2}"
    @argcheck !Ferrite.isclosed(ch) "The ConstraintHandler is already closed"
    displacementDofs(p) = [fieldOffset(dh, f) + _nComponents(dh, f) * (p - 1) + c for f in fields for c in 1:_nComponents(dh, f)]
    @argcheck length(displacementDofs(1)) == 3 "pass one 3D vector field or three scalar displacement fields"

    edge, inward = _edgeAndInwardControlPoints(dh, side, first(fields))
    base_ip = Ferrite.get_base_interpolation(_fieldInterpolation(dh, first(fields)))
    grevilleU, grevilleV = (toVector(TinyGismo.greville(TinyGismo.knots(_maybeDeref(base_ip.basis), d))) for d in 1:2)
    greville = [[u, v] for v in grevilleV for u in grevilleU]
    geometry = Ferrite.get_grid(dh).geometry

    # Only collected here: a control point at a corner gets a condition from each side, and
    # those have to be solved together in `close!`.
    conditions = get!(Dict{Vector{Int}, Vector{Vector{Pair{Int, Float64}}}}, _normalConditions, ch)
    for (e, p) in zip(edge, inward)
        d₀, d₁ = displacementDofs(e), displacementDofs(p)
        J = toMatrix(TinyGismo.jacobian(geometry, greville[e]))
        n = normalize(cross(J[:, 1], J[:, 2]))
        push!(get!(Vector{Vector{Pair{Int, Float64}}}, conditions, d₁), [d₁ .=> n; d₀ .=> -n])
    end
    return nothing
end

# Pending conditions of the tuple form, keyed by the dofs of the control point next to the side.
# Kept outside the ConstraintHandler since Ferrite's struct has no room for them.
const _normalConditions = WeakKeyDict{ConstraintHandler, Dict{Vector{Int}, Vector{Vector{Pair{Int, Float64}}}}}()

function Ferrite.close!(ch::ConstraintHandler{<:IGADofHandler})
    conditions = get(_normalConditions, ch, nothing)
    if !isnothing(conditions)
        # Sorted, so the result does not depend on the order of the fixEdgeTangent! calls
        for d₁ in sort!(collect(keys(conditions)))
            _addNormalConstraints!(ch, sort!(conditions[d₁]))
        end
        delete!(_normalConditions, ch)
    end
    return invoke(Ferrite.close!, Tuple{ConstraintHandler}, ch)
end

#=
Ferrite picks the constrained dof of each condition and untangles them in `close!`, but rejects
conditions without a free dof and linearly dependent ones. Both occur at corners: the control
point next to a side can lie on a held edge, and two sides with the same normal give the same
condition, which then only holds if both edges move alike along n.
=#
function _addNormalConstraints!(ch::ConstraintHandler, conditions)
    freeDofs = unique(first(e) for c in conditions for e in c if !haskey(ch.dofmapping, first(e)))
    rows = zeros(0, length(freeDofs))
    for condition in conditions
        row = zeros(1, length(freeDofs))
        for (dof, coefficient) in condition
            j = findfirst(==(dof), freeDofs)
            isnothing(j) || (row[j] += coefficient)
        end
        candidate = [rows; row]
        rank(candidate; atol = sqrt(eps())) > size(rows, 1) || continue
        rows = candidate
        add!(ch, AffineConstraint(condition, 0.0))
    end
    return ch
end

# The control points on `side` and the row next to it. The first parametric direction runs
# fastest, so that row is one control point away for left/right and one whole row for bottom/top.
function _edgeAndInwardControlPoints(dh::IGADofHandler, side, field_name::Symbol)
    base_ip = Ferrite.get_base_interpolation(_fieldInterpolation(dh, field_name))
    @argcheck !(_maybeDeref(base_ip.basis) isa TinyGismo.HierarchicalBasis) "fixEdgeTangent! does not support hierarchical (THB) patches"
    edge = collect(edgeControlPoints(dh, side, field_name))
    side_index = isa(side, Symbol) ? _side_idx_for_symbol(side) : side
    stride = side_index <= 2 ? 1 : Int(TinyGismo.size(_maybeDeref(base_ip.basis))) ÷ length(edge)
    return edge, edge .+ (isodd(side_index) ? stride : -stride)
end
