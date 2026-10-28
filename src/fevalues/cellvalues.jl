# The active basis functions depend on the knot span, so the reference-cell quadrature
# points are remapped into the span's parameter space and the interpolation is re-evaluated
# there. Everything after that (mapping, `detJdV`, hessians incl. embedded elements) is
# Ferrite's generic `reinit!`.
function _remap_points(points, span::KnotSpanWrapper{dim}) where {dim}
    return [Vec{dim}(Tuple(ref_to_param(ξ, span))) for ξ in points]
end

function _precompute_values!(cv_or_fv, points)
    foreach(fv -> Ferrite.precompute_values!(fv, points), Ferrite.get_fun_values(cv_or_fv))
    Ferrite.precompute_values!(Ferrite.get_geo_mapping(cv_or_fv), points)
    return nothing
end

function Ferrite.reinit!(cv::CellValues, cc::CellCache{<:Any, <:IGAGrid})
    span = getspan(cc.grid, Ferrite.cellid(cc))
    _precompute_values!(cv, _remap_points(Ferrite.getpoints(Ferrite.get_quadrature_rule(cv)), span))

    Ferrite.reinit!(cv, nothing, cc.coords)

    # Volume mapping from the reference cell to the knot span
    rV = areaOfKnotSpan(span) / 2^length(span.lower)
    detJdV = Ferrite.getdetJdVs(cv)
    detJdV === nothing || (detJdV .*= rV)
    return nothing
end
