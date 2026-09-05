# ==============================================================================
# Boundary facet sets
# ==============================================================================

"""
    getfacetset(grid::IGAGrid, name::String)

Return the set of boundary facets (as `Ferrite.FacetIndex` entries) of `grid` on the
parametric boundary `name`. Since an `IGAGrid` is a structured tensor-product grid of knot
spans, the facet sets do not need to be stored; they are computed on demand from the
knot-span corners.

Valid names are `"left"`, `"right"` (and, for two-dimensional grids, `"bottom"`, `"top"`),
referring to the sides of the parametric domain. The returned set can be used with
`FacetIterator` and `FacetValues` in the usual Ferrite way:

# Example
```julia
for facet in FacetIterator(dh, getfacetset(grid, "right"))
    reinit!(facetvalues, facet)
    # ...
end
```
"""
function Ferrite.getfacetset(grid::IGAGrid{sdim, 2}, name::String) where {sdim}
    name in ("left", "right", "bottom", "top") || throw(
        ArgumentError("Unknown facet set \"$name\"; expected \"left\", \"right\", \"bottom\" or \"top\"")
    )
    lo = grid.knotSpans[1].lower
    hi = grid.knotSpans[end].upper
    set = Ferrite.OrderedSet{FacetIndex}()
    for (i, kv) in enumerate(grid.knotSpans)
        # Ferrite facet numbering of `RefQuadrilateral`: 1=bottom, 2=right, 3=top, 4=left
        if name == "bottom" && kv.lower[2] == lo[2]
            push!(set, FacetIndex(i, 1))
        elseif name == "right" && kv.upper[1] == hi[1]
            push!(set, FacetIndex(i, 2))
        elseif name == "top" && kv.upper[2] == hi[2]
            push!(set, FacetIndex(i, 3))
        elseif name == "left" && kv.lower[1] == lo[1]
            push!(set, FacetIndex(i, 4))
        end
    end
    return set
end

function Ferrite.getfacetset(grid::IGAGrid{sdim, 1}, name::String) where {sdim}
    name in ("left", "right") || throw(
        ArgumentError("Unknown facet set \"$name\"; expected \"left\" or \"right\"")
    )
    facet = name == "left" ? FacetIndex(1, 1) : FacetIndex(getncells(grid), 2)
    return Ferrite.OrderedSet{FacetIndex}((facet,))
end

# Ferrite restricts the convenience constructor to `Union{Grid, AbstractDofHandler}`, so
# the plain-grid variant needs to be provided explicitly for `IGAGrid`.
function Ferrite.FacetIterator(
        grid::IGAGrid, set::Ferrite.AbstractVecOrSet{FacetIndex},
        flags::UpdateFlags = UpdateFlags()
    )
    return Ferrite.FacetIterator(Ferrite.FacetCache(grid, flags), set)
end

# ==============================================================================
# FacetValues reinit! for IGA interpolations
# ==============================================================================

# Scaling of the reference-facet quadrature weights to the parametric facet measure,
# i.e. the 1D analogue of `areaOfKnotSpan` for the facet's tangential direction(s).
function facetScaleOfKnotSpan(knotSpan::KnotSpanWrapper{2}, facet_nr::Int)
    l, u = knotSpan.lower, knotSpan.upper
    # Facets 1 (bottom) and 3 (top) run along parametric direction 1, facets 2 (right)
    # and 4 (left) along direction 2.
    return facet_nr in (1, 3) ? (u[1] - l[1]) / 2 : (u[2] - l[2]) / 2
end

# In 1D the facets are points, so no weight scaling is needed.
facetScaleOfKnotSpan(::KnotSpanWrapper{1}, ::Int) = 1.0

function Ferrite.reinit!(
        fv::FacetValues{FV}, cc::CellCache, facet_nr::Int
    ) where {
        FV <:
        Ferrite.FunctionValues{Order, IP} where {
            Order,
            IP <: Union{IGAInterpolation, VectorizedInterpolation{<:Any, <:Any, <:Any, <:IGAInterpolation}},
        },
    }
    Ferrite.set_current_facet!(fv, facet_nr)
    fun_values = Ferrite.get_fun_values(fv)
    geo_mapping = Ferrite.get_geo_mapping(fv)

    knotSpan = cc.grid.knotSpans[Ferrite.cellid(cc)]
    dim = Ferrite.getrefdim(fun_values.ip)

    # Remap the facet quadrature points into the active knot span, as in the cell reinit!
    points = [Vec{dim}(Tuple(ref_to_param(ξ, knotSpan))) for ξ in Ferrite.getpoints(fv.fqr, facet_nr)]
    Ferrite.precompute_values!(fun_values, points)
    Ferrite.precompute_values!(geo_mapping, points)

    current_coefs = getfield.(cc.grid.nodes[cc.nodes], :x)

    rF = facetScaleOfKnotSpan(knotSpan, facet_nr)

    @inbounds for (q_point, w) in pairs(Ferrite.getweights(fv.fqr, facet_nr))
        mapping = Ferrite.calculate_mapping(geo_mapping, q_point, current_coefs)
        J = Ferrite.getjacobian(mapping)
        # `J` is the Jacobian w.r.t. the parametric (knot) coordinates; the affine map from
        # the reference facet to the parametric facet is accounted for by `rF`.
        weight_norm = Ferrite.weighted_normal(J, Ferrite.getrefshape(geo_mapping.ip), facet_nr)
        detJ = norm(weight_norm)
        detJ > 0.0 || Ferrite.throw_detJ_not_pos(detJ)
        fv.detJdV[q_point] = detJ * w * rF
        fv.normals[q_point] = weight_norm / detJ
        Ferrite.apply_mapping!(fun_values, q_point, mapping, cc)
    end
    return nothing
end
