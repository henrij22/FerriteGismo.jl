"""
    IGADofHandler(grid::IGAGrid)

A degree-of-freedom handler for [`IGAGrid`](@ref)s. It mirrors the interface of the
regular `Ferrite.DofHandler` — fields are added with `add!` and the handler is finalized
with `close!` — but distributes the degrees of freedom according to the globally numbered
spline control points instead of Lagrange nodes.

Internally it wraps a `Ferrite.DofHandler` and forwards most property accesses to it, so
that downstream Ferrite functionality (`allocate_matrix`, `CellIterator`, `celldofs!`,
`ndofs`, …) keeps working.

Its `SubDofHandler`s are owned by the `IGADofHandler` itself rather than by the wrapped
`Ferrite.DofHandler`, so that `sdh.dh` refers back to the IGA handler. That is what lets
`CellCache`/`CellIterator` over a `SubDofHandler` pick the IGA specializations, and it is
required for subdomain-wise iteration.

# Example
```julia
dh = IGADofHandler(grid)
add!(dh, :u, ip)
close!(dh)
```

!!! warning
    The `IGADofHandler` relies on internal Ferrite functionality and only supports a
    single subdofhandler / cell type. See also [`fixEdge!`](@ref) for prescribing
    homogeneous Dirichlet conditions on the parametric edges.
"""
mutable struct IGADofHandler{dim, G <: IGAGrid{dim}} <: Ferrite.AbstractDofHandler
    dh::DofHandler{dim, G}
    # Owned (not forwarded): these are bound to the IGADofHandler, see the docstring above
    subdofhandlers::Vector{SubDofHandler{IGADofHandler{dim, G}}}
    field_offsets::Vector{Int}
end

const _IGA_OWN_FIELDS = (:dh, :subdofhandlers, :field_offsets)

# Access to wrapped dofhandler components
function Base.getproperty(dh::IGADofHandler, sym::Symbol)
    return (sym in _IGA_OWN_FIELDS) ? getfield(dh, sym) : getproperty(getfield(dh, :dh), sym)
end
function Base.setproperty!(dh::IGADofHandler, sym::Symbol, v)
    return (sym in _IGA_OWN_FIELDS) ? setfield!(dh, sym, v) : setproperty!(getfield(dh, :dh), sym, v)
end

function IGADofHandler(grid::G) where {dim, G <: IGAGrid{dim}}
    sdhs = SubDofHandler{IGADofHandler{dim, G}}[]
    return IGADofHandler{dim, G}(DofHandler(grid), sdhs, Int[])
end

# Mirrors `Ferrite.ndofs_per_cell(::DofHandler)` over the owned `SubDofHandler`s; forwarding
# to the wrapped handler would hit its (now empty) list.
function Ferrite.ndofs_per_cell(dh::IGADofHandler)
    if length(dh.subdofhandlers) > 1
        error("There are more than one subdofhandler. Use `ndofs_per_cell(dh, cellid::Int)` instead.")
    end
    @assert length(dh.subdofhandlers) != 0
    return @inbounds ndofs_per_cell(dh.subdofhandlers[1])
end
function Ferrite.ndofs_per_cell(dh::IGADofHandler, cell::Int)
    sdhidx = dh.cell_to_subdofhandler[cell]
    sdhidx ∉ 1:length(dh.subdofhandlers) && return 0 # only defined on a subdomain
    return ndofs_per_cell(dh.subdofhandlers[sdhidx])
end

#=
Mirrors `Ferrite.add!(::DofHandler, ...)` but creates the `SubDofHandler` on the
`IGADofHandler` itself. Delegating to the wrapped handler (as this used to do) would bind
the `SubDofHandler` to the plain `Ferrite.DofHandler`, so `sdh.dh` would hand back an
unwrapped handler and the IGA `CellCache`/`celldofs!` specializations would be bypassed.
=#
function Ferrite.add!(dh::IGADofHandler, name::Symbol, ip::Interpolation)
    @assert !Ferrite.isclosed(dh)
    celltype = getcelltype(Ferrite.get_grid(dh))
    @assert isconcretetype(celltype)
    if isempty(dh.subdofhandlers)
        # Create a new SubDofHandler for all cells
        sdh = SubDofHandler(dh, 1:getncells(Ferrite.get_grid(dh)))
    elseif length(dh.subdofhandlers) == 1
        # Add to existing SubDofHandler (if it covers all cells)
        sdh = dh.subdofhandlers[1]
        if length(sdh.cellset) != getncells(Ferrite.get_grid(dh))
            error("can not add field to IGADofHandler with a SubDofHandler for a subregion")
        end
    else
        error("can not add field to IGADofHandler with multiple SubDofHandlers")
    end
    add!(sdh, name, ip)
    return dh
end

# `find_field`/`getfieldinterpolation` iterate `dh.subdofhandlers`, which now live on the
# IGA handler — the wrapped `Ferrite.DofHandler` has none.
function Ferrite.find_field(dh::IGADofHandler, field_name::Symbol)
    for (sdh_idx, sdh) in pairs(dh.subdofhandlers)
        field_idx = Ferrite._find_field(sdh, field_name)
        !isnothing(field_idx) && return (sdh_idx, field_idx)
    end
    return error("Did not find field :$field_name in IGADofHandler (existing fields: $(Ferrite.getfieldnames(dh))).")
end

function Ferrite.getfieldinterpolation(dh::IGADofHandler, field_idxs::NTuple{2, Int})
    sdh_idx, field_idx = field_idxs
    return dh.subdofhandlers[sdh_idx].field_interpolations[field_idx]
end

function Ferrite.close!(dh::IGADofHandler)
    @assert !Ferrite.isclosed(dh)

    # Collect the global field names
    empty!(dh.field_names)
    for sdh in dh.subdofhandlers, name in sdh.field_names
        name in dh.field_names || push!(dh.field_names, name)
    end

    # One global dof block per field, shared by every SubDofHandler. The hierarchical
    # subdomains partition the *elements* of one patch, so a control point active in two of
    # them must come out as the same dof in both.
    empty!(dh.field_offsets)
    fieldOffsets = Dict{Symbol, Int}()
    nextdof = 0
    for name in dh.field_names
        fieldOffsets[name] = nextdof
        push!(dh.field_offsets, nextdof)
        nextdof += _fieldDofCount(dh, name)
    end

    ncells = getncells(dh.grid)
    cellDofs = [Int[] for _ in 1:ncells]

    for (sdhi, sdh) in pairs(dh.subdofhandlers)
        _collect_cell_dofs_iga!(cellDofs, dh, sdh, fieldOffsets)

        sdh.ndofs_per_cell = 0
        for ip in sdh.field_interpolations
            sdh.ndofs_per_cell += Ferrite.get_base_interpolation(ip).nbasefuns * Ferrite.n_components(ip)
        end

        for ci in sdh.cellset
            @assert dh.cell_to_subdofhandler[ci] == 0
            dh.cell_to_subdofhandler[ci] = sdhi
        end
    end

    # Flatten in cell order, recording where each cell starts. The offsets cannot be a fixed
    # stride: subdomains may carry different numbers of dofs per cell.
    empty!(dh.cell_dofs)
    for i in 1:ncells
        dh.cell_dofs_offset[i] = length(dh.cell_dofs) + 1
        append!(dh.cell_dofs, cellDofs[i])
    end

    dh.ndofs = maximum(dh.cell_dofs; init = 0)
    dh.closed = true

    return dh
end

# Field lookup across subdomains. Which interpolation describes a field, and where its dof
# block starts, are properties of the field and its basis rather than of any one subdomain,
# so `only(dh.subdofhandlers)` is too strict for a hierarchically grouped patch.

"""
    _globalFieldIndex(dh::IGADofHandler, name::Symbol) -> Int

Index of `name` in `dh.field_names`, which is also its index into `dh.field_offsets`.
"""
function _globalFieldIndex(dh::IGADofHandler, name::Symbol)
    idx = findfirst(==(name), dh.field_names)
    idx === nothing && error("Did not find field :$name in IGADofHandler (existing fields: $(Ferrite.getfieldnames(dh))).")
    return idx
end

"""
    fieldOffset(dh::IGADofHandler, field_name::Symbol) -> Int

The dof before `field_name`'s global block, i.e. `u[fieldOffset(dh, f) + ncomp * (cp - 1) + c]`.

Each field occupies one contiguous block, shared by every subdomain carrying it. Use this
rather than indexing `dh.field_offsets` with a `SubDofHandler`-local field index; the two
coincide only for a single subdomain.
"""
fieldOffset(dh::IGADofHandler, field_name::Symbol) = dh.field_offsets[_globalFieldIndex(dh, field_name)]

"""
    _fieldInterpolation(dh::IGADofHandler, name::Symbol) -> Interpolation

The interpolation describing field `name`, taken from whichever subdomain defines it.
"""
function _fieldInterpolation(dh::IGADofHandler, name::Symbol)
    for sdh in dh.subdofhandlers
        idx = Ferrite._find_field(sdh, name)
        idx === nothing || return sdh.field_interpolations[idx]
    end
    return error("Did not find field :$name in IGADofHandler (existing fields: $(Ferrite.getfieldnames(dh))).")
end

# Total dofs of one field. Subdomains share a field's global block, so they must agree on its
# size; numbering it from the first alone would silently corrupt the rest.
function _fieldDofCount(dh::IGADofHandler, name::Symbol)
    count = nothing
    for sdh in dh.subdofhandlers
        idx = Ferrite._find_field(sdh, name)
        idx === nothing && continue
        ip = sdh.field_interpolations[idx]
        n = Int(TinyGismo.size(Ferrite.get_base_interpolation(ip).basis)) * Ferrite.n_components(ip)
        if count === nothing
            count = n
        elseif count != n
            error(
                "field :$name is described by different interpolations on different " *
                    "SubDofHandlers ($count vs $n dofs). Subdomains of one patch share a " *
                    "field's dofs, so they must use the same basis."
            )
        end
    end
    count === nothing && error("field :$name is not defined on any SubDofHandler")
    return count
end

function _collect_cell_dofs_iga!(
        cellDofs::Vector{Vector{Int}}, dh::IGADofHandler, sdh::SubDofHandler,
        fieldOffsets::Dict{Symbol, Int}
    )
    grid = dh.grid
    actives = gsMatrix{Int32}()

    for (name, ip) in zip(sdh.field_names, sdh.field_interpolations)
        offset = fieldOffsets[name]
        basis = Ferrite.get_base_interpolation(ip).basis
        ncomp = Ferrite.n_components(ip)

        for ci in sdh.cellset
            # Read from the grid rather than re-derived from the basis.
            activeInCell = _activeIn(basis, grid.knotSpans[ci], actives)

            if ncomp == 1
                append!(cellDofs[ci], activeInCell .+ offset)
            else
                # Interweave dofs u1x, u1y, u1z, u2x, ...
                for a in activeInCell, c in 1:ncomp
                    push!(cellDofs[ci], offset + ncomp * (a - 1) + c)
                end
            end
        end
    end
    return cellDofs
end

#=
Ferrite's `allocate_matrix`/`init_sparsity_pattern`/`add_cell_entries!` are all typed on the
concrete `Ferrite.DofHandler` and read `dh.subdofhandlers`, so they cannot be reached by
forwarding to the wrapped handler (which owns none) nor by dispatching on the IGA handler.
The pattern is therefore assembled here directly from the cell dofs.
=#
function Ferrite.allocate_matrix(dh::IGADofHandler)
    @assert Ferrite.isclosed(dh)
    n = ndofs(dh)
    nnz_per_row = 2 * maximum(Ferrite.ndofs_per_cell, dh.subdofhandlers; init = 1)
    sp = Ferrite.SparsityPattern(n, n; nnz_per_row = nnz_per_row)
    for sdh in dh.subdofhandlers
        cellDofs = zeros(Int, Ferrite.ndofs_per_cell(sdh))
        for cellid in sdh.cellset
            Ferrite.celldofs!(cellDofs, dh, cellid)
            for col in cellDofs, row in cellDofs
                Ferrite.add_entry!(sp, row, col)
            end
        end
    end
    return Ferrite.allocate_matrix(sp)
end
Ferrite.get_grid(dh::IGADofHandler) = Ferrite.get_grid(dh.dh)
Ferrite.getfieldnames(dh::IGADofHandler) = Ferrite.getfieldnames(dh.dh)
