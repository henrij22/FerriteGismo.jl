# ==============================================================================
# L2 projection onto the spline space
#
# Quadrature point data (stresses, …) is discontinuous and has to be projected onto the
# spline space before it can be visualized. This is a self-contained implementation on
# top of the public Ferrite API — `CellValues`, `allocate_matrix`, `start_assemble` /
# `assemble!` — rather than a specialization of Ferrite's `L2Projector` internals.
# ==============================================================================

"""
    L2ProjectorIGA <: Ferrite.AbstractProjector

IGA counterpart of Ferrite's `L2Projector`, used to project (typically discontinuous)
quadrature-point data such as stresses onto a continuous spline field for postprocessing
and VTK export.

It is constructed through `Ferrite.L2Projector(ip::IGAInterpolation, grid::IGAGrid)` and
stores the Cholesky factorization of the spline mass matrix together with an internal
[`IGADofHandler`](@ref) carrying one scalar field, and the interpolations used to build it.
Use `project` to obtain the projected control-point values, `write_projection` to export
them and `evaluate_at_grid_nodes` to read them off the export mesh.
"""
struct L2ProjectorIGA{IP <: IGAInterpolation, GIP <: IGAInterpolation, DH <: IGADofHandler, F} <: Ferrite.AbstractProjector
    M_cholesky::F
    dh::DH
    ip::IP
    ip_geo::GIP
    qr_lhs::QuadratureRule
end

# Quadrature sufficient for integrating the mass matrix of a spline basis of degree `order`
_massQuadratureRule(::IGAInterpolation{shape, order}) where {shape, order} = QuadratureRule{shape}(order + 1)
_massQuadratureRule(ip::VectorizedInterpolation{<:Any, <:Any, <:Any, <:IGAInterpolation}) = _massQuadratureRule(ip.ip)

# The geometry of a patch is described by its own basis, so the geometric interpolation of
# the values is always the grid's, independently of the space projected onto.
function _geometryInterpolation(grid::IGAGrid, ip::Interpolation)
    return IGAInterpolation{getrefshape(ip)}(TinyGismo.basis(grid.geometry))
end

"""
    L2Projector(ip::IGAInterpolation, grid::IGAGrid; qr_lhs = <order + 1>, set = all cells)

Create a projector that projects quadrature-point data onto the spline space `ip` over
`grid`, see [`L2ProjectorIGA`](@ref).

A vector-valued interpolation (`ip^dim`) may be passed as well; the projection always
happens component-wise onto the scalar base interpolation. `set` is accepted for
compatibility with Ferrite's signature, but an `IGAGrid` is a single patch, so it has to
cover all cells.
"""
function Ferrite.L2Projector(
        ip::IP,
        grid::GRID;
        qr_lhs::QuadratureRule = _massQuadratureRule(ip),
        set = 1:getncells(grid),
        geom_ip = nothing
    ) where {
        GRID <: IGAGrid,
        IP <: Union{IGAInterpolation, VectorizedInterpolation{<:Any, <:Any, <:Any, <:IGAInterpolation}},
    }
    geom_ip === nothing ||
        @warn("Providing geom_ip is deprecated, the geometric interpolation of the cells will always be used")
    length(set) == getncells(grid) ||
        error("An IGAGrid is a single patch: the L2Projector has to cover all of its cells")

    ip_fun = ip isa VectorizedInterpolation ? ip.ip : ip
    ip_geo = _geometryInterpolation(grid, ip_fun)

    dh = IGADofHandler(grid)
    add!(dh, :_, ip_fun)
    close!(dh)

    M = _assembleMassMatrix(dh, ip_fun, ip_geo, qr_lhs)
    return L2ProjectorIGA(cholesky(Symmetric(M)), dh, ip_fun, ip_geo, qr_lhs)
end

function _assembleMassMatrix(dh::IGADofHandler, ip, ip_geo, qr::QuadratureRule)
    cv = CellValues(qr, ip, ip_geo; update_gradients = false)
    M = allocate_matrix(dh)
    assembler = start_assemble(M)
    n = getnbasefunctions(ip)
    Me = zeros(n, n)

    for cell in CellIterator(dh)
        reinit!(cv, cell)
        fill!(Me, 0.0)
        for q_point in 1:getnquadpoints(cv)
            dΩ = getdetJdV(cv, q_point)
            for i in 1:n
                Nᵢ = shape_value(cv, q_point, i)
                for j in 1:n
                    Me[i, j] += Nᵢ * shape_value(cv, q_point, j) * dΩ
                end
            end
        end
        assemble!(assembler, celldofs(cell), Me)
    end
    return M
end

"""
    project(proj::L2ProjectorIGA, vars, qr_rhs::QuadratureRule)

Project the quadrature-point data `vars` onto the spline space of `proj` and return the
resulting control-point values, one per basis function of the projection space.

`vars` is indexed by **cell id** and holds, per cell, one value per quadrature point of
`qr_rhs`, in the order in which the rule visits them. Scalars and tensors are supported;
tensor components are projected independently.
"""
function Ferrite.project(proj::L2ProjectorIGA, vars::AbstractVector{TC}, qr_rhs::QuadratureRule) where {T, TC <: AbstractVector{T}}
    length(vars) == getncells(Ferrite.get_grid(proj.dh)) ||
        error("vars is indexed by the cellid: length(vars) != number of cells")

    M = T <: Tensors.AbstractTensor ? Tensors.n_components(Tensors.get_base(T)) : 1
    f = _assembleProjectionRHS(proj, vars, qr_rhs, M)

    projected = proj.M_cholesky \ f
    makeT(vals) = T <: Tensors.AbstractTensor ? T(Tuple(vals)) : vals[1]
    return T[makeT(row) for row in eachrow(projected)]
end

# f = ∫ N ⋅ x̂ dΩ, one column per component of the (possibly tensorial) data
function _assembleProjectionRHS(proj::L2ProjectorIGA, vars, qr_rhs::QuadratureRule, M::Int)
    cv = CellValues(qr_rhs, proj.ip, proj.ip_geo; update_gradients = false)
    n = getnbasefunctions(proj.ip)
    f = zeros(ndofs(proj.dh), M)

    component(x::Tensors.AbstractTensor, i) = x.data[i]
    component(x::Number, _) = x

    for cell in CellIterator(proj.dh)
        cellVars = vars[cellid(cell)]
        length(cellVars) == getnquadpoints(qr_rhs) ||
            error("The number of variables per cell doesn't match the number of quadrature points")
        reinit!(cv, cell)
        dofs = celldofs(cell)

        for q_point in 1:getnquadpoints(cv)
            dΩ = getdetJdV(cv, q_point)
            qpVars = cellVars[q_point]
            for i in 1:n
                Nᵢ = shape_value(cv, q_point, i)
                for j in 1:M
                    f[dofs[i], j] += Nᵢ * component(qpVars, j) * dΩ
                end
            end
        end
    end
    return f
end

"""
    evaluate_at_grid_nodes(proj::L2ProjectorIGA, vals::AbstractVector; subdivision::Int = 1)

Evaluate projected values `vals` (the output of `project`) at the nodes of the export mesh,
see [`exportGrid`](@ref).
"""
function Ferrite.evaluate_at_grid_nodes(proj::L2ProjectorIGA, vals::AbstractVector; subdivision::Int = 1)
    return evaluateAtExportNodes(proj.dh, vals, :_; subdivision)
end
