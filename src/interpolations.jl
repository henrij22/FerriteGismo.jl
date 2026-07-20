"""
    IGAInterpolation{shape}(basis::gsBasis)

A `Ferrite.ScalarInterpolation` that wraps a G+Smo spline `basis` (B-spline or NURBS)
so that it can be used like any other Ferrite interpolation, e.g. when constructing
`CellValues` or when adding a field to an [`IGADofHandler`](@ref).

`shape` is the Ferrite reference shape of the elements, e.g. `RefQuadrilateral` for a
two-dimensional tensor-product basis or `RefLine` in one dimension. The polynomial
`order` is deduced from the degree of `basis`.

Because the active basis functions differ from knot span to knot span, an
`IGAInterpolation` is *mutable* and keeps track of the currently active element; this is
updated automatically during `reinit!` and should not be modified by user code.

Vector-valued fields are created in the usual Ferrite way, e.g. `ip^2` for a
two-component field.

# Example
```julia
ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
qr = QuadratureRule{RefQuadrilateral}(2)
cv = CellValues(qr, ip, ip)
```
"""
mutable struct IGAInterpolation{shape, order, B, dim} <: ScalarInterpolation{shape, order}
    basis::B
    nbasefuns::Int
    currentElement::KnotSpanWrapper{dim}
end

struct IGAMapping end
Ferrite.mapping_type(::IGAInterpolation) = IGAMapping()

function IGAInterpolation{shape}(basis::BB) where {shape <: Ferrite.AbstractRefShape, BB}
    dim = Ferrite.getrefdim(shape)

    currentElement = KnotSpanWrapper{dim}(first(knotSpans(basis)))
    if dim == 1
        order = TinyGismo.degree(basis)
        nbasefuns = numActive(basis) # might be a problem with tensorbasis
    else
        order = maximum(ntuple(i -> TinyGismo.degree(basis, i), dim))
        out = gsMatrix{Int32}()
        active!(basis, Vector(currentElement.center), out)
        nbasefuns = TinyGismo.rows(out)
    end

    return IGAInterpolation{shape, order, BB, dim}(
        basis, nbasefuns, currentElement
    )
end

Ferrite.getnbasefunctions(ip::IGAInterpolation) = ip.nbasefuns

function Ferrite.reference_shape_values!(
        values::AbstractMatrix, ip::IP, qr_points::AbstractVector{<:Vec{rdim}}
    ) where {rdim, IP <: IGAInterpolation}
    @boundscheck checkbounds(values, 1:getnbasefunctions(ip))

    valsRaw = gsMatrix()
    for (qp, ξref) in pairs(qr_points)
        ξ = ref_to_param(ξref, ip.currentElement.lower, ip.currentElement.upper)
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
    for (qp, ξref) in pairs(qr_points)
        ξ = ref_to_param(ξref, ip.currentElement.lower, ip.currentElement.upper)
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
    for (qp, ξref) in pairs(qr_points)
        ξ = ref_to_param(ξref, ip.currentElement.lower, ip.currentElement.upper)
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


# THis depends actually on the cell
Ferrite.conformity(::IGAInterpolation) = Ferrite.H1Conformity()

# ==============================================================================
# Vector-valued (vectorized) IGA interpolations
# ==============================================================================
# A `VectorizedInterpolation` (e.g. `ip^2`) wrapping an `IGAInterpolation` is used for
# vector fields such as displacements. Ferrite's generic fallback evaluates such
# interpolations one basis function at a time via automatic differentiation, which requires
# a scalar `reference_shape_value(::IGAInterpolation, ξ, i)` method that does not exist for
# spline bases (they are evaluated for all active functions at once). We therefore provide
# specialized batch methods that evaluate the scalar spline basis once and expand the result
# to the vectorized layout (node-major, component-minor: u1x, u1y, u2x, u2y, ...).

# The vectorized interpolation must use the IGA mapping, not the identity mapping.
Ferrite.mapping_type(::VectorizedInterpolation{<:Any, <:Any, <:Any, <:IGAInterpolation}) = IGAMapping()

_iga_base(ip::IGAInterpolation) = ip
_iga_base(ip::VectorizedInterpolation{<:Any, <:Any, <:Any, <:IGAInterpolation}) = ip.ip

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
