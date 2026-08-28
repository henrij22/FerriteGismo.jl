"""
    KnotSpanWrapper{dim}

Internal, immutable description of a single knot span (Bézier element) of a tensor-product
spline basis. It caches the `center`, `lower` and `upper` corner of the knot span in
parameter space as `Vec{dim}` values, converted once from the underlying
`TinyGismo.KnotSpan`.

`KnotSpanWrapper`s back the cells of an [`IGAGrid`](@ref) and are used during `reinit!` to
remap quadrature points to parameter space (via `FerriteGismo.ref_to_param`) before an
[`IGAInterpolation`](@ref) evaluates them, and to compute knot-span areas.
"""
struct KnotSpanWrapper{dim}
    center::Vec{dim}
    lower::Vec{dim}
    upper::Vec{dim}
    function KnotSpanWrapper{dim}(knotSpan::TinyGismo.KnotSpan) where {dim}
        return new{dim}(toVec(Vec{dim}, centerPoint(knotSpan)), toVec(Vec{dim}, lowerCorner(knotSpan)), toVec(Vec{dim}, upperCorner(knotSpan)))
    end

    # Hierarchical bases have no `TinyGismo.KnotSpan`s to wrap -- their elements come out of
    # `elementBoxes` as plain corner coordinates, so build the wrapper from those instead.
    function KnotSpanWrapper{dim}(lower::AbstractVector, upper::AbstractVector) where {dim}
        lo = Vec{dim}(ntuple(i -> Float64(lower[i]), dim))
        up = Vec{dim}(ntuple(i -> Float64(upper[i]), dim))
        return new{dim}((lo + up) / 2, lo, up)
    end
end

toVec(::Type{Vec{dim}}, vec::gsVector{T}) where {dim, T} = Vec{dim, T}(toVector(vec))
# ==============================================================================
# Utility functions
# ==============================================================================

extractCoefs(coefs::Matrix) = extractCoefs(coefs, typeof(zero(Vec{2})))

function extractCoefs(coefs::Matrix, ::Type{T}) where {T <: AbstractVector}
    return [coefs[i, :] |> T for i in 1:size(coefs, 1)]
end

"""
    interpolate(basis::gsBasis, u::AbstractVector, x; offset = 0)

Evaluate the spline field defined by the coefficient/dof vector `u` over the G+Smo `basis`
at the parametric point `x` (either a `Ferrite.Vec` or a plain vector).

The active basis functions at `x` are looked up and combined with the corresponding
entries of `u`. `offset` is added to the (1-based) coefficient indices before indexing into
`u`, which is used internally to address a particular field inside a global solution
vector. The element type of the result follows `eltype(u)`, so both scalar and
tensor-valued fields are supported.
"""
function interpolate(basis::gsBasis, u::AbstractVector, x::Vec; offset = 0)
    return interpolate(basis, u, Vector(x); offset)
end

function interpolate(basis::gsBasis, u::AbstractVector{T}, x::AbstractVector; offset = 0) where {T}
    result = zero(T)
    shape_values = gsMatrix()
    actives = gsMatrix{Int32}()

    active!(basis, x, actives)
    TinyGismo.eval!(basis, x, shape_values)

    shape_values = toVector(shape_values)
    actives = toVector(actives)
    actives .+= offset

    for i in eachindex(shape_values)
        result += shape_values[i] * u[actives[i]]
    end
    return result
end

_maybeDeref(obj::TinyGismo.CxxWrap.CxxWrapCore.ConstCxxRef) = obj[]
_maybeDeref(obj) = obj

"""
    _elementSpans(basis, ::Val{dim}) -> Vector{KnotSpanWrapper{dim}}

The elements of `basis` in parameter space, one wrapper per element.

This is the seam over how a basis enumerates its elements. Tensor-product bases go through
`TinyGismo.knotSpans`; hierarchical bases cannot, because G+Smo's hierarchical domain
iterator is not safely copyable and TinyGismo therefore refuses `knotSpans` for them. They
use `TinyGismo.elementBoxes` instead, which returns the same information as owned data.
"""
_elementSpans(basis, v::Val) = _elementSpans(_maybeDeref(basis), v)

_elementSpans(basis::gsBasis, ::Val{dim}) where {dim} = map(KnotSpanWrapper{dim}, knotSpans(basis))

function _elementSpans(basis::TinyGismo.HierarchicalBasis, ::Val{dim}) where {dim}
    boxes = toMatrix(elementBoxes(basis))
    return [
        KnotSpanWrapper{dim}(view(boxes, 1:dim, j), view(boxes, (dim + 1):(2dim), j))
            for j in axes(boxes, 2)
    ]
end

"""
    _activeIn(basis, span::KnotSpanWrapper, out::gsMatrix{Int32}) -> Vector{Int32}

The basis functions active on one element, as global 1-based indices.

Evaluated at the element *center* rather than a corner: on a hierarchical mesh a corner can
sit on the boundary between two levels, where the active set is not the one belonging to
this element.
"""
function _activeIn(basis, span::KnotSpanWrapper, out = gsMatrix{Int32}())
    active!(basis, Vector(span.center), out)
    return toVector(out)
end
