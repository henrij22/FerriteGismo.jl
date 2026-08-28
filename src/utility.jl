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
