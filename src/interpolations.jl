"""
    IGAInterpolation{shape}(basis::gsBasis)

A `Ferrite.ScalarInterpolation` that wraps a G+Smo spline `basis` (B-spline or NURBS)
so that it can be used like any other Ferrite interpolation, e.g. when constructing
`CellValues` or when adding a field to an [`IGADofHandler`](@ref).

`shape` is the Ferrite reference shape of the elements, e.g. `RefQuadrilateral` for a
two-dimensional tensor-product basis or `RefLine` in one dimension. The polynomial
`order` is deduced from the degree of `basis`.

Since the active basis functions depend on the current knot span, `reinit!` remaps
quadrature points into that span's parameter space (see
[`FerriteGismo.KnotSpanWrapper`](@ref)) before evaluating `ip`. `ip` itself holds no
per-cell state and is as safe to share as any other Ferrite interpolation.

Vector-valued fields are created in the usual Ferrite way, e.g. `ip^2` for a
two-component field.

# Example
```julia
ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
qr = QuadratureRule{RefQuadrilateral}(2)
cv = CellValues(qr, ip, ip)
```
"""
struct IGAInterpolation{shape, order, B} <: ScalarInterpolation{shape, order}
    basis::B
    nbasefuns::Int
end

struct IGAMapping end
Ferrite.mapping_type(::IGAInterpolation) = IGAMapping()

function IGAInterpolation{shape}(basis::BB) where {shape <: Ferrite.AbstractRefShape, BB}
    dim = Ferrite.getrefdim(shape)

    if dim == 1
        order = TinyGismo.degree(basis)
        nbasefuns = numActive(basis) # might be a problem with tensorbasis
    else
        order = maximum(ntuple(i -> TinyGismo.degree(basis, i), dim))
        # nbasefuns is the same for every knot span, so any one will do as a probe point
        firstSpan = KnotSpanWrapper{dim}(first(knotSpans(basis)))
        out = gsMatrix{Int32}()
        active!(basis, Vector(firstSpan.center), out)
        nbasefuns = TinyGismo.rows(out)
    end

    return IGAInterpolation{shape, order, BB}(basis, nbasefuns)
end

Ferrite.getnbasefunctions(ip::IGAInterpolation) = ip.nbasefuns

# `qr_points` are already remapped into parameter space by `reinit!`, not the canonical
# [-1, 1]^d reference-cell points the name suggests.
function Ferrite.reference_shape_values!(
        values::AbstractMatrix, ip::IP, qr_points::AbstractVector{<:Vec{rdim}}
    ) where {rdim, IP <: IGAInterpolation}
    @boundscheck checkbounds(values, 1:getnbasefunctions(ip))

    valsRaw = gsMatrix()
    for (qp, ξ) in pairs(qr_points)
        eval!(ip.basis, Vector(ξ), valsRaw)
        for i in 1:getnbasefunctions(ip)
            values[i, qp] = valsRaw[i, 1]
        end
    end
    return
end

function Ferrite.reference_shape_gradients_and_values!(
        gradients::AbstractMatrix, values::AbstractMatrix, ip::IP, qr_points::AbstractVector{<:Vec{rdim}}
    ) where {rdim, IP <: IGAInterpolation}
    @boundscheck checkbounds(gradients, 1:getnbasefunctions(ip))
    @boundscheck checkbounds(values, 1:getnbasefunctions(ip))

    valsRaw = gsMatrix()
    derivsRaw = gsMatrix()
    for (qp, ξ) in pairs(qr_points)
        eval!(ip.basis, Vector(ξ), valsRaw)
        deriv!(ip.basis, Vector(ξ), derivsRaw)

        @inbounds for i in 1:getnbasefunctions(ip)
            values[i, qp] = valsRaw[i, 1]
            gradients[i, qp] = Vec{rdim}(j -> (derivsRaw[i * (rdim) - (rdim - j), 1]))
        end
    end
    return
end

function Ferrite.reference_shape_hessians_gradients_and_values!(
        hessians::AbstractMatrix,
        gradients::AbstractMatrix, values::AbstractMatrix, ip::IP, qr_points::AbstractVector{<:Vec{rdim}}
    ) where {rdim, IP <: IGAInterpolation}
    @boundscheck checkbounds(gradients, 1:getnbasefunctions(ip))
    @boundscheck checkbounds(values, 1:getnbasefunctions(ip))

    valsRaw = gsMatrix()
    derivsRaw = gsMatrix()
    derivs2Raw = gsMatrix()
    for (qp, ξ) in pairs(qr_points)
        eval!(ip.basis, Vector(ξ), valsRaw)
        deriv!(ip.basis, Vector(ξ), derivsRaw)
        deriv2!(ip.basis, Vector(ξ), derivs2Raw)

        @inbounds for i in 1:getnbasefunctions(ip)
            values[i, qp] = valsRaw[i, 1]
            gradients[i, qp] = Vec{rdim}(j -> (derivsRaw[i * (rdim) - (rdim - j), 1]))
            hessians[i, qp] = SymmetricTensor{2, rdim}(
                rdim == 2 ? (derivs2Raw[3i - 2, 1], derivs2Raw[3i, 1], derivs2Raw[3i - 1, 1]) : (derivs2Raw[i, 1],)
            )
        end
    end
    return
end


# This depends actually on the cell
Ferrite.conformity(::IGAInterpolation) = Ferrite.H1Conformity()

# ==============================================================================
# Vector-valued (vectorized) IGA interpolations
# ==============================================================================
# Vector fields (e.g. `ip^2`) can't use Ferrite's per-basis-function AD fallback, since that
# needs a scalar `reference_shape_value(::IGAInterpolation, ξ, i)` which spline bases lack
# (all active functions are evaluated at once). We instead evaluate the scalar basis once and
# expand it to the vectorized layout (node-major, component-minor: u1x, u1y, u2x, u2y, ...).

# The vectorized interpolation must use the IGA mapping, not the identity mapping.
Ferrite.mapping_type(::VectorizedInterpolation{<:Any, <:Any, <:Any, <:IGAInterpolation}) = IGAMapping()

function Ferrite.reference_shape_values!(
        values::AbstractMatrix, ipv::VectorizedInterpolation{vdim, shape, order, IP},
        qr_points::AbstractVector{<:Vec}
    ) where {vdim, shape, order, IP <: IGAInterpolation}
    ip = ipv.ip
    T = eltype(eltype(values))
    n = getnbasefunctions(ip)
    scalar_vals = Matrix{T}(undef, n, length(qr_points))
    Ferrite.reference_shape_values!(scalar_vals, ip, qr_points)

    for qp in eachindex(qr_points)
        for I in 1:getnbasefunctions(ipv)
            i0, c0 = divrem(I - 1, vdim)
            i, c = i0 + 1, c0 + 1
            v = scalar_vals[i, qp]
            values[I, qp] = Vec{vdim, T}(j -> j == c ? v : zero(T))
        end
    end
    return
end

function Ferrite.reference_shape_gradients_and_values!(
        gradients::AbstractMatrix, values::AbstractMatrix,
        ipv::VectorizedInterpolation{vdim, shape, order, IP},
        qr_points::AbstractVector{<:Vec{rdim}}
    ) where {vdim, rdim, shape, order, IP <: IGAInterpolation}
    ip = ipv.ip
    T = eltype(eltype(values))
    n = getnbasefunctions(ip)
    scalar_vals = Matrix{T}(undef, n, length(qr_points))
    scalar_grads = Matrix{Vec{rdim, T}}(undef, n, length(qr_points))
    Ferrite.reference_shape_gradients_and_values!(scalar_grads, scalar_vals, ip, qr_points)

    for qp in eachindex(qr_points)
        for I in 1:getnbasefunctions(ipv)
            i0, c0 = divrem(I - 1, vdim)
            i, c = i0 + 1, c0 + 1
            v = scalar_vals[i, qp]
            g = scalar_grads[i, qp]
            values[I, qp] = Vec{vdim, T}(j -> j == c ? v : zero(T))
            gradients[I, qp] = Tensor{2, vdim, T}((a, b) -> a == c ? g[b] : zero(T))
        end
    end
    return
end

"""
    interpolate(ip::IGAInterpolation, u::AbstractVector, x; offset = 0)
    interpolate(ip::VectorizedInterpolation{vdim, …, <:IGAInterpolation}, u::AbstractVector, x; offset = 0)

Evaluate the field described by the dof vector `u` and the interpolation `ip` at the
parametric point `x`.

The scalar version simply forwards to the `gsBasis` method above. For a vector-valued field
(e.g. `ip^2`) the interleaved dof numbering of the [`IGADofHandler`](@ref) (`u1x, u1y, u2x,
u2y, …`) is taken into account and a `Vec{vdim}` is returned. `offset` is the offset of the
field inside the global dof vector.
"""
interpolate(ip::IGAInterpolation, u::AbstractVector, x; offset = 0) = interpolate(_maybeDeref(ip.basis), u, x; offset)

function interpolate(
        ip::VectorizedInterpolation{vdim, shape, order, <:IGAInterpolation},
        u::AbstractVector{T}, x::Vec; offset = 0
    ) where {vdim, shape, order, T}
    return interpolate(ip, u, Vector(x); offset)
end

function interpolate(
        ip::VectorizedInterpolation{vdim, shape, order, <:IGAInterpolation},
        u::AbstractVector{T}, x::AbstractVector; offset = 0
    ) where {vdim, shape, order, T}
    basis = _maybeDeref(ip.ip.basis)

    shape_values_gs = gsMatrix()
    actives_gs = gsMatrix{Int32}()
    active!(basis, x, actives_gs)
    TinyGismo.eval!(basis, x, shape_values_gs)
    shape_values = toVector(shape_values_gs)
    actives = toVector(actives_gs)

    result = zero(Vec{vdim, T})
    for k in eachindex(shape_values)
        p = actives[k]
        result += shape_values[k] * Vec{vdim, T}(ntuple(c -> u[offset + vdim * (p - 1) + c], vdim))
    end
    return result
end
