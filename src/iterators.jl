#=
Ferrite 1.7 changed `CellCache` from a mutable struct with `cellid::Int` to an immutable one
whose `cellid` is a one-element array (so that caches can be built on views/GPU arrays).
The constructor below is shared with older versions by feeding it whichever the loaded
Ferrite expects; `Ferrite.cellid(cc)` is the version-agnostic way to read it back.
=#
const CELLCACHE_HAS_ARRAY_CELLID = pkgversion(Ferrite) >= v"1.7"
_initial_cellid() = CELLCACHE_HAS_ARRAY_CELLID ? [-1] : -1

function Ferrite.CellCache(grid::IGAGrid{sdim, rdim}, flags::UpdateFlags = UpdateFlags()) where {sdim, rdim}
    N = Ferrite.nnodes_per_cell(grid, 1) # nodes and coords will be resized in `reinit!`
    nodes = zeros(Int, N)
    coords = zeros(Vec{sdim, Ferrite.get_coordinate_eltype(grid)}, N)
    return CellCache(flags, grid, _initial_cellid(), nodes, coords, nothing, Int[])
end

function Ferrite.CellCache(
        dh::DH, flags::UpdateFlags = UpdateFlags()
    ) where {DH <: IGADofHandler{dim, G}} where {dim, G}
    n = ndofs_per_cell(dh.subdofhandlers[1]) # dofs and coords will be resized in `reinit!`
    N = Ferrite.nnodes_per_cell(Ferrite.get_grid(dh), 1)
    nodes = zeros(Int, N)
    coords = zeros(Vec{dim, Ferrite.get_coordinate_eltype(Ferrite.get_grid(dh))}, N)
    celldofs = zeros(Int, n)
    return CellCache(flags, Ferrite.get_grid(dh), _initial_cellid(), nodes, coords, dh, celldofs)
end

function Ferrite.celldofs!(global_dofs::AbstractVector{Int}, dh::IGADofHandler, i::Int)
    @assert Ferrite.isclosed(dh)
    @assert length(global_dofs) == Ferrite.ndofs_per_cell(dh, i)
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

#=
Iterating a `SubDofHandler` has to go through the IGA cache too. Without these, Ferrite's
generic `CellCache(sdh)`/`CellIterator(sdh)` would be used and the IGA `reinit!`/`celldofs!`
path bypassed. These are what make subdomain-wise iteration work for IGA, which is required
once several `SubDofHandler`s share one handler.
=#
function Ferrite.CellCache(sdh::SubDofHandler{<:IGADofHandler}, flags::UpdateFlags = UpdateFlags())
    return CellCache(sdh.dh, flags)
end

function Ferrite.CellIterator(sdh::SubDofHandler{<:IGADofHandler}, flags::UpdateFlags = UpdateFlags())
    return Ferrite.CellIterator(CellCache(sdh, flags), sdh.cellset)
end
