# The basis depends on the knot span, so the quadrature points are remapped and the
# interpolation re-evaluated before Ferrite's generic `reinit!`.
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

    rV = areaOfKnotSpan(span) / 2^length(span.lower)
    detJdV = Ferrite.getdetJdVs(cv)
    detJdV === nothing || (detJdV .*= rV)
    return nothing
end
