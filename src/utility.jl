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
