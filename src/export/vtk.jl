Ferrite.cell_to_vtkcell(::Type{IGACell{1}}) = Ferrite.VTKCellTypes.VTK_LINE
Ferrite.cell_to_vtkcell(::Type{IGACell{2}}) = Ferrite.VTKCellTypes.VTK_QUAD
Ferrite.cell_to_vtkcell(::Type{IGACell{3}}) = Ferrite.VTKCellTypes.VTK_HEXAHEDRON

function Ferrite.create_vtk_griddata(grid::IGAGrid{sdim}) where {sdim}
    uniformGrid = parameterSpaceGrid(grid)
    coords, cls = Ferrite.create_vtk_griddata(uniformGrid)
    coordsPhysical = zeros(sdim, size(coords, 2))
    for i in axes(coords, 2)
        coordsPhysical[:, i] = toPhysical(grid, coords[:, i])
    end
    return coordsPhysical, cls
end

function Ferrite.VTKGridFile(filename::String, dh::IGADofHandler; kwargs...)
    return VTKGridFile(filename, dh.dh; kwargs...)
end

function Ferrite._assemble_L2_matrix(dh::IGADofHandler, qrs_lhs::Vector{<:QuadratureRule})
    M = Symmetric(allocate_matrix(dh))
    assembler = start_assemble(M)
    for (sdh, qr_lhs) in zip(dh.subdofhandlers, qrs_lhs)
        ip_fun = only(sdh.field_interpolations)
        ip_geo = Ferrite.default_geometric_interpolation(ip_fun)
        cv = CellValues(qr_lhs, ip_fun, ip_geo; update_gradients = false)
        Ferrite._assemble_L2_matrix!(assembler, cv, sdh)
    end
    return M
end

mutable struct L2ProjectorIGA <: Ferrite.AbstractProjector
    M_cholesky::Any #::SuiteSparse.CHOLMOD.Factor{Float64}
    dh::IGADofHandler
    qrs_lhs::Vector{<:QuadratureRule}
    qrs_rhs::Vector{<:QuadratureRule}
end

function L2ProjectorIGA(grid::IGAGrid)
    dh = IGADofHandler(grid)
    return L2ProjectorIGA(nothing, dh, QuadratureRule[], QuadratureRule[])
end

function Ferrite.L2Projector(
        ip::IGAInterpolation,
        grid::IGAGrid;
        qr_lhs::QuadratureRule = _mass_qr(ip),
        set = Ferrite.OrderedSet(1:getncells(grid)),
        geom_ip = nothing
    )
    geom_ip === nothing ||
        @warn("Providing geom_ip is deprecated, the geometric interpolation of the cells with always be used")
    proj = L2ProjectorIGA(grid)
    add!(proj, set, ip; qr_lhs, qr_rhs = nothing)
    close!(proj)
    return proj
end
Ferrite.isclosed(proj::L2ProjectorIGA) = Ferrite.isclosed(proj.dh.dh)

function Ferrite.add!(
        proj::L2ProjectorIGA, set::Ferrite.AbstractVecOrSet{Int}, ip::Interpolation;
        qr_rhs::Union{QuadratureRule, Nothing}, qr_lhs::QuadratureRule = _mass_qr(ip)
    )
    # Validate user input
    Ferrite.isclosed(proj) && error("The L2Projector is already closed")
    if qr_rhs !== nothing
        Ferrite.getrefshape(ip) == Ferrite.getrefshape(qr_rhs) ||
            error("The reference shape of the interpolation and the qr_rhs must be the same")
        push!(proj.qrs_rhs, qr_rhs)
    end
    Ferrite.getrefshape(ip) == Ferrite.getrefshape(qr_lhs) ||
        error("The reference shape of the interpolation and the qr_lhs must be the same")

    # sdh = SubDofHandler(proj.dh, set)
    add!(proj.dh, :_, ip isa VectorizedInterpolation ? ip.ip : ip)
    push!(proj.qrs_lhs, qr_lhs)

    return proj
end

"""
    close!(proj::L2Projector)

Close `proj` which assembles and calculates the left-hand-side of the projection equation, before doing a Cholesky factorization
of the mass-matrix.
"""
function Ferrite.close!(proj::L2ProjectorIGA)
    close!(proj.dh)
    M = Ferrite._assemble_L2_matrix(proj.dh, proj.qrs_lhs)
    proj.M_cholesky = cholesky(Symmetric(M))
    return proj
end

function Ferrite.project(
        p::L2ProjectorIGA, vars::Union{AbstractVector, AbstractDict}, qr_rhs::QuadratureRule
    )
    length(p.dh.subdofhandlers) == 1 ||
        error("For multiple domains, provide the right-hand-side quadrature rule to the L2Projector")
    return Ferrite._project(p, vars, [qr_rhs])
end

function Ferrite._project(
        proj::L2ProjectorIGA,
        vars::Union{AbstractVector{TC}, AbstractDict{Int, TC}},
        qrs_rhs::Vector{<:QuadratureRule}
    ) where {
        T <: Union{Number, AbstractTensor}, TC <: AbstractVector{T},
    }

    # Sanity checks for user input
    Ferrite.isclosed(proj) || error("The L2Projector is not closed")
    length(qrs_rhs) == 0 &&
        error("The right-hand-side quadrature rule must be provided, unless already given to the L2Projector")
    length(qrs_rhs) == length(proj.dh.subdofhandlers) ||
        error("Number of qrs_rhs must match the number of `add!`ed sets")
    for (qr_rhs, sdh) in zip(qrs_rhs, proj.dh.subdofhandlers)
        if getrefshape(qr_rhs) !== getrefshape(getcelltype(sdh))
            error("Reference shape of quadrature rule and cells doesn't match. Please ensure that `qrs_rhs` has the same order as sets are added to the L2Projector")
        end
    end
    # Catch if old input-style giving vars indexed by the set index, instead of the cell id
    if isa(vars, AbstractVector) && length(vars) != getncells(Ferrite.get_grid(proj.dh))
        error("vars is indexed by the cellid, not the index in the set: length(vars) != number of cells")
    end

    M = T <: AbstractTensor ? Tensors.n_components(Tensors.get_base(T)) : 1

    return Ferrite._project(proj, qrs_rhs, vars, M, T)::Vector{T}
end

function Ferrite._project(
        proj::L2ProjectorIGA, qrs_rhs::Vector{<:QuadratureRule},
        vars::Union{AbstractVector, AbstractDict}, M::Integer, ::Type{T}
    ) where {T}
    f = zeros(ndofs(proj.dh), M)
    for (sdh, qr_rhs) in zip(proj.dh.subdofhandlers, qrs_rhs)
        ip_fun = only(sdh.field_interpolations)
        ip_geo = Ferrite.default_geometric_interpolation(ip_fun)
        cv = CellValues(qr_rhs, ip_fun, ip_geo; update_gradients = false)
        Ferrite.assemble_proj_rhs!(f, cv, sdh, vars)
    end

    # solve for the projected nodal values
    projected_vals = proj.M_cholesky \ f

    # Recast to original input type
    make_T(vals) = T <: AbstractTensor ? T(Tuple(vals)) : vals[1]
    return T[make_T(x) for x in Base.eachrow(projected_vals)]
end

# function Ferrite.assemble_proj_rhs!(
#         f::Matrix, cellvalues::CellValues{FV}, sdh::SubDofHandler, vars::Union{AbstractVector, AbstractDict}) where {FV <: Ferrite.FunctionValues{Order, IP} where {Order, IP <: IGAInterpolation}}
#     # Assemble the multi-column rhs, f = ∭( v ⋅ x̂ )dΩ
#     # The number of columns corresponds to the length of the data-tuple in the tensor x̂.
#     M = size(f, 2)
#     n = getnbasefunctions(cellvalues)
#     fe = zeros(n, M)
#     nqp = getnquadpoints(cellvalues)

#     get_data(x::AbstractTensor, i) = x.data[i]
#     get_data(x::Number, _) = x

#     dh = IGADofHandler(sdh.dh, Int[])

#     ## Assemble contributions from each cell
#     for cell in CellIterator(dh)
#         fill!(fe, 0)
#         cell_vars = vars[cellid(cell)]
#         length(cell_vars) == nqp ||
#             error("The number of variables per cell doesn't match the number of quadrature points")
#         reinit!(cellvalues, cell)

#         for q_point in 1:nqp
#             dΩ = getdetJdV(cellvalues, q_point)
#             qp_vars = cell_vars[q_point]
#             for i in 1:n
#                 v = shape_value(cellvalues, q_point, i)
#                 for j in 1:M
#                     fe[i, j] += v * get_data(qp_vars, j) * dΩ
#                 end
#             end
#         end

#         # Assemble cell contribution
#         for (num, dof) in enumerate(celldofs(cell))
#             f[dof, :] += fe[num, :]
#         end
#     end
#     return
# end

function Ferrite.write_projection(vtk::VTKGridFile, proj::L2ProjectorIGA, vals, name)
    if Ferrite.write_discontinuous(vtk)
        data = Ferrite.evaluate_at_discontinuous_vtkgrid_nodes(
            proj.dh, vals, only(getfieldnames(proj.dh)), vtk.cellnodes
        )
    else
        data = Ferrite._evaluate_at_grid_nodes(proj, vals, Val(true))::Matrix #=vtk=#
    end
    Ferrite._vtk_write_node_data(
        vtk.vtk, data, name; component_names = Ferrite.component_names(eltype(vals))
    )
    return vtk
end

function Ferrite.evaluate_at_grid_nodes(proj::L2ProjectorIGA, vals::AbstractVector)
    return Ferrite._evaluate_at_grid_nodes(proj, vals, Val(false))
end

# Numbers can be handled by the method for DofHandler
function Ferrite._evaluate_at_grid_nodes(proj::L2ProjectorIGA, vals::AbstractVector{<:Number}, vtk)
    return Ferrite._evaluate_at_grid_nodes(proj.dh, vals, only(Ferrite.getfieldnames(proj.dh)), vtk)
end

function Ferrite._evaluate_at_grid_nodes(
        proj::L2ProjectorIGA, vals::AbstractVector{S},
        ::Val{vtk}
    ) where {order, dim, T, M, S <: Union{Tensor{order, dim, T, M}, SymmetricTensor{order, dim, T, M}}, vtk}
    dh = proj.dh
    uniformGrid = parameterSpaceGrid(dh.grid)
    # The internal dofhandler in the projector is a scalar field, but the values in vals
    # can be any tensor field, however, the number of dofs should always match the length of vals
    @assert ndofs(dh) == length(vals)
    if vtk
        nout = S <: Vec{2} ? 3 : M # Pad 2D Vec to 3D
        data = fill(T(NaN), nout, Ferrite.getnnodes(uniformGrid))
    else
        data = fill(T(NaN) * zero(S), Ferrite.getnnodes(Ferrite.get_grid(dh)))
    end
    for sdh in dh.subdofhandlers
        ip = Ferrite.getfieldinterpolation(sdh, 1)
        _evaluate_at_grid_nodes_iga!(data, sdh, vals, ip, dh.field_offsets[1])
    end
    return data
end
