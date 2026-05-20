function _side_idx_for_symbol(side::Symbol)::Int
    if side == :left
        return 1
    elseif side == :right
        return 2
    elseif side == :bottom
        return 3
    elseif side == :top
        return 4
    else
        throw(ArgumentError("Unknown side $side"))
    end
end


function fixEdge!(dh::IGADofHandler, ch::ConstraintHandler, side::Union{Symbol, Integer}, field_name::Symbol)
    return fixEdge!(dh, ch, side, Ferrite.find_field(only(dh.dh.subdofhandlers), field_name))
end

function fixEdge!(dh::IGADofHandler, ch::ConstraintHandler, side::Union{Symbol, Integer}, field_index::Int)
    side_index = isa(side, Symbol) ? _side_idx_for_symbol(side) : side
    field_ip = Ferrite.getfieldinterpolation(dh.dh, (1, field_index))
    indices = toMatrix(boundary(_maybeDeref(field_ip.basis), side_index)) .+ dh.field_offsets[field_index]

    for i in indices
        Ferrite.add_prescribed_dof!(ch, i, 0.0)
    end
    return nothing
end
