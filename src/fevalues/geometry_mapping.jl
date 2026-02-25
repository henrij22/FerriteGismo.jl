Ferrite.required_geo_diff_order(::IGAMapping, fun_diff_order::Int) = fun_diff_order

# Identity mapping
@inline function Ferrite.apply_mapping!(
        funvals::Ferrite.FunctionValues{0}, ::IGAMapping, q_point::Int, mapping_values, args...
    )
    @inbounds for j in 1:getnbasefunctions(funvals)
        funvals.Nx[j, q_point] = funvals.Nξ[j, q_point]
    end
    return nothing
end

@inline function Ferrite.apply_mapping!(
        funvals::Ferrite.FunctionValues{1}, ::IGAMapping, q_point::Int, mapping_values, args...
    )
    Jinv = Ferrite.calculate_Jinv(Ferrite.getjacobian(mapping_values))
    @inbounds for j in 1:getnbasefunctions(funvals)
        funvals.Nx[j, q_point] = funvals.Nξ[j, q_point]
        # funvals.dNdx[j, q_point] = funvals.dNdξ[j, q_point] ⋅ Jinv # TODO via Tensors.jl
        funvals.dNdx[j, q_point] = Ferrite.dothelper(funvals.dNdξ[j, q_point], Jinv)
    end
    return nothing
end

@inline function Ferrite.apply_mapping!(
        funvals::Ferrite.FunctionValues{2}, ::IGAMapping, q_point::Int, mapping_values, args...
    )
    Jinv = Ferrite.calculate_Jinv(Ferrite.getjacobian(mapping_values))
    sdim, rdim = size(Jinv)

    H = rdim == sdim ? Ferrite.gethessian(mapping_values) : nothing
    is_vector_valued = first(funvals.Nx) isa Vec
    Jinv_otimesu_Jinv = is_vector_valued && difforder > 1 ? Ferrite.otimesu(Jinv, Jinv) : nothing
    @inbounds for j in 1:getnbasefunctions(funvals)
        funvals.Nx[j, q_point] = funvals.Nξ[j, q_point]
        dNdx = Ferrite.dothelper(funvals.dNdξ[j, q_point], Jinv)
        funvals.dNdx[j, q_point] = dNdx

        if rdim == sdim
            if is_vector_valued
                d2Ndx2 = (funvals.d2Ndξ2[j, q_point] - dNdx ⋅ H) ⊡ Jinv_otimesu_Jinv
            else
                d2Ndx2 = Jinv' ⋅ (funvals.d2Ndξ2[j, q_point] - dNdx ⋅ H) ⋅ Jinv
            end
            funvals.d2Ndx2[j, q_point] = d2Ndx2
        end
    end
    return nothing
end

@inline function Ferrite.calculate_mapping(
        geo_mapping::Ferrite.GeometryMapping{2, IP}, q_point::Int,
        x::AbstractVector{<:Vec}
    ) where {IP <: IGAInterpolation}
    J = zero(Ferrite.otimes_returntype(eltype(x), eltype(geo_mapping.dMdξ)))
    sdim, rdim = size(J)
    # (rdim != sdim) && error("hessian for embedded elements not implemented (rdim=$rdim, sdim=$sdim)")
    H = rdim == sdim ? zero(Ferrite.otimes_returntype(eltype(x), eltype(geo_mapping.d2Mdξ2))) : nothing
    @inbounds for j in 1:Ferrite.getngeobasefunctions(geo_mapping)
        J += Ferrite.otimes_helper(x[j], geo_mapping.dMdξ[j, q_point])
        if sdim == rdim
            H += x[j] ⊗ geo_mapping.d2Mdξ2[j, q_point]
        end
    end
    return Ferrite.MappingValues(J, H)
end
