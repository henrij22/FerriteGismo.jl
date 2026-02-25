function Ferrite.CellCache(grid::IGAGrid{sdim, rdim}, flags::UpdateFlags = UpdateFlags()) where {sdim, rdim}
    N = Ferrite.nnodes_per_cell(grid, 1) # nodes and coords will be resized in `reinit!`
    nodes = zeros(Int, N)
    coords = zeros(Vec{sdim, Ferrite.get_coordinate_eltype(grid)}, N)
    return CellCache(flags, grid, -1, nodes, coords, nothing, Int[])
end

function Ferrite.CellCache(
        dh::DH, flags::UpdateFlags = UpdateFlags()
    ) where {DH <: IGADofHandler{dim, G}} where {dim, G}
    n = ndofs_per_cell(dh.subdofhandlers[1]) # dofs and coords will be resized in `reinit!`
    N = Ferrite.nnodes_per_cell(Ferrite.get_grid(dh), 1)
    nodes = zeros(Int, N)
    coords = zeros(Vec{dim, Ferrite.get_coordinate_eltype(Ferrite.get_grid(dh))}, N)
    celldofs = zeros(Int, n)
    return CellCache(flags, Ferrite.get_grid(dh), -1, nodes, coords, dh, celldofs)
end

function Ferrite.celldofs!(global_dofs, dh::IGADofHandler, i::Int)
    @assert Ferrite.isclosed(dh)
    @assert length(global_dofs) == Ferrite.ndofs_per_cell(dh)
    unsafe_copyto!(global_dofs, 1, dh.cell_dofs, dh.cell_dofs_offset[i], length(global_dofs))
    return global_dofs
end

function Ferrite.CellIterator(
        gridordh::Union{IGAGrid, IGADofHandler},
        set::Union{Ferrite.IntegerCollection, Nothing} = nothing,
        flags::UpdateFlags = UpdateFlags()
    )
    if set === nothing
        grid = gridordh isa IGADofHandler ? Ferrite.get_grid(gridordh) : gridordh
        set = 1:getncells(grid)
    end
    if gridordh isa IGADofHandler
        # TODO: Since the CellCache is resizeable this is not really necessary to check
        #       here, but might be useful to catch slow code paths?
        Ferrite._check_same_celltype(Ferrite.get_grid(gridordh), set)
    end
    return Ferrite.CellIterator(CellCache(gridordh, flags), set)
end
function Ferrite.CellIterator(gridordh::Union{IGAGrid, IGADofHandler}, flags::UpdateFlags)
    return CellIterator(gridordh, nothing, flags)
end
