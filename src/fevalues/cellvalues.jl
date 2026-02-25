function Ferrite.reinit!(
        cv::CellValues{FV}, cell::CellCache
    ) where {
        FV <:
        Ferrite.FunctionValues{Order, IP} where {Order, IP <: IGAInterpolation},
    }
    cv.fun_values.ip.currentElement = cell.grid.knotSpans[cell.cellid]
    cv.geo_mapping.ip.currentElement = cell.grid.knotSpans[cell.cellid]

    Ferrite.precompute_values!(cv.fun_values, Ferrite.getpoints(cv.qr))
    Ferrite.precompute_values!(cv.geo_mapping, Ferrite.getpoints(cv.qr))

    geo_mapping = cv.geo_mapping
    fun_values = cv.fun_values

    # Ferrite.check_reinit_sdim_consistency(:IGACellValues, Ferrite.shape_gradient_type(cv), eltype(x))
    if cell === nothing && reinit_needs_cell(cv)
        throw(ArgumentError("The cell::AbstractCell input is required to reinit! non-identity function mappings"))
    end

    current_coefs = getfield.(cell.grid.nodes[cell.nodes], :x)
    dim = Ferrite.getrefdim(cv.fun_values.ip)

    # Volume Mapping from parent to knot span
    rV = areaOfKnotSpan(geo_mapping.ip.currentElement) / (2^dim)

    @inbounds for (q_point, w) in enumerate(Ferrite.getweights(cv.qr))
        mapping = Ferrite.calculate_mapping(geo_mapping, q_point, current_coefs)
        Ferrite._update_detJdV!(cv.detJdV, q_point, w, mapping)
        @inbounds cv.detJdV[q_point] *= rV
        Ferrite.apply_mapping!(fun_values, q_point, mapping, cell)
    end
    return nothing
end
