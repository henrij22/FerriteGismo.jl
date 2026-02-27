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

# TODO
# function Ferrite.default_geometric_interpolation(
#         ip::IGAInterpolation{shape}
#     ) where {shape <: Ferrite.AbstractRefShape{dim}}
#     return VectorizedInterpolation{dim}(
#         IGAInterpolation{shape, order}(
#             Gismo.basis(ip.geometry), ip.geometry
#         )
#     )
# end
# function Ferrite.default_geometric_interpolation(
#         ::IGAInterpolation{
#             vdim, shape, order,
#             IGAInterpolation{shape, order},
#         }
#     ) where {vdim, order, dim, shape <: Ferrite.AbstractRefShape{dim}}
#     return VectorizedInterpolation{dim}(
#         IGAInterpolation{shape, order}(
#             Gismo.basis(ip.geometry), ip.geometry
#         )
#     )
# end

# THis depends actually on the cell
Ferrite.conformity(::IGAInterpolation) = Ferrite.H1Conformity()
