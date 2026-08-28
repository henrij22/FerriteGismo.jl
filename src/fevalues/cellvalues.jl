function Ferrite.reinit!(
        cv::CellValues{FV}, cell::CellCache
    ) where {
        FV <:
        Ferrite.FunctionValues{Order, IP} where {
            Order,
            IP <: Union{IGAInterpolation, VectorizedInterpolation{<:Any, <:Any, <:Any, <:IGAInterpolation}},
        },
    }
    # `IGAInterpolation` carries no per-cell state: instead of switching some "active
    # element" stored on the interpolation, the reference-cell quadrature points are
    # remapped into the active knot span's parameter-space rectangle here, before ever
    # being handed to `ip`. This is what makes `ip` (and, transitively, `cv`'s per-task
    # scratch buffers aside, `cv` itself) safe to share across cells and tasks.
    knotSpan = cell.grid.knotSpans[cell.cellid]
    points = map(ξ -> ref_to_param(ξ, knotSpan), Ferrite.getpoints(cv.qr))

    Ferrite.precompute_values!(cv.fun_values, points)
    Ferrite.precompute_values!(cv.geo_mapping, points)

    geo_mapping = cv.geo_mapping
    fun_values = cv.fun_values

    # Ferrite.check_reinit_sdim_consistency(:IGACellValues, Ferrite.shape_gradient_type(cv), eltype(x))
    if cell === nothing && reinit_needs_cell(cv)
        throw(ArgumentError("The cell::AbstractCell input is required to reinit! non-identity function mappings"))
    end

    current_coefs = getfield.(cell.grid.nodes[cell.nodes], :x)
    dim = Ferrite.getrefdim(cv.fun_values.ip)

    # Volume Mapping from parent to knot span
    rV = areaOfKnotSpan(knotSpan) / (2^dim)

    @inbounds for (q_point, w) in enumerate(Ferrite.getweights(cv.qr))
        mapping = Ferrite.calculate_mapping(geo_mapping, q_point, current_coefs)
        Ferrite._update_detJdV!(cv.detJdV, q_point, w, mapping)
        @inbounds cv.detJdV[q_point] *= rV
        Ferrite.apply_mapping!(fun_values, q_point, mapping, cell)
    end
    return nothing
end
