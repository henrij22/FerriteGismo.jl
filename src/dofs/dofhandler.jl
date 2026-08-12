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

    # Set initial values
    nextdof = 1  # next free dof to distribute

    for (sdhi, sdh) in pairs(dh.subdofhandlers)
        nextdof = _close_subdofhandler_iga!(dh, sdh, sdhi, nextdof)
    end
    dh.ndofs = maximum(dh.cell_dofs; init = 0)
    dh.closed = true

    return dh
end

function _close_subdofhandler_iga!(
        dh::IGADofHandler, sdh::SubDofHandler, sdh_index::Int,
        nextdof::Int
    )
    dof_offsets = Int[]
    sdh.ndofs_per_cell = 0

    v = Vector{Int}[]
    for i in 1:getncells(dh.grid)
        push!(v, Int[])
    end

    current_dof_counter = 0
    actives = gsMatrix{Int32}()
    for ip in sdh.field_interpolations
        kvs = map(FerriteGismo.KnotSpanWrapper{get_rdim(dh.dh.grid)}, knotSpans(Ferrite.get_base_interpolation(ip).basis))
        ncomp = Ferrite.n_components(ip)
        basisSize = TinyGismo.size(Ferrite.get_base_interpolation(ip).basis) * ncomp

        for (i, kv) in enumerate(kvs)
            if ip isa VectorizedInterpolation
                active!(Ferrite.get_base_interpolation(ip).basis, Vector(kv.lower), actives)
                activeInCell = toVector(actives)

                # Interweave dofs u1x, u1y, u1z, u2x, ...
                # `current_dof_counter` is a dof offset, so it is added outside the
                # component interleaving (matching fixEdge! and the grid-node evaluation).
                for j in 1:length(activeInCell)
                    for c in 1:ncomp
                        push!(v[i], current_dof_counter + ncomp * (activeInCell[j] - 1) + c)
                    end
                end
            else
                active!(ip.basis, Vector(kv.lower), actives)
                activeInCell = toVector(actives) .+ current_dof_counter
                push!(v[i], activeInCell...)
            end
        end
        push!(dh.field_offsets, current_dof_counter)
        current_dof_counter += basisSize
    end

    for i in eachindex(dh.grid.knotSpans)
        push!(dh.cell_dofs, v[i]...)
    end

    # Check: Is this redundant information?
    for ip in sdh.field_interpolations
        sdh.ndofs_per_cell += Ferrite.get_base_interpolation(ip).nbasefuns * Ferrite.n_components(ip)
    end

    dh.cell_dofs_offset .= collect(1:(sdh.ndofs_per_cell):(getncells(dh.grid) * sdh.ndofs_per_cell))

    for ci in sdh.cellset
        @assert dh.cell_to_subdofhandler[ci] == 0
        dh.cell_to_subdofhandler[ci] = sdh_index
    end

    return maximum(dh.cell_dofs)
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

function Ferrite._evaluate_at_grid_nodes(
        dh::IGADofHandler{sdim}, u::AbstractVector{T}, fieldname::Symbol,
        ::Val{vtk} = Val(false)
    ) where {T, vtk, sdim}
    uniformGrid = parameterSpaceGrid(dh.grid)

    # Make sure the field exists
    fieldname ∈ Ferrite.getfieldnames(dh) || error("Field $fieldname not found.")
    # Figure out the return type (scalar or vector)
    field_idx = Ferrite.find_field(dh, fieldname)
    ip = Ferrite.getfieldinterpolation(dh, field_idx)
    if vtk
        # VTK output of solution field (or L2 projected scalar data)
        n_c = Ferrite.n_components(ip)
        vtk_dim = n_c == 2 ? 3 : n_c # VTK wants vectors padded to 3D
        # Float32 is the smallest float type supported by VTK
        TT = promote_type(T, Float32)
        data = fill!(Matrix{TT}(undef, vtk_dim, Ferrite.getnnodes(uniformGrid)), NaN)
    else
        # Just evaluation at grid nodes
        RT = typeof(Ferrite.function_value_init(ip, u))
        data = fill(T(NaN) * zero(RT), getnnodes(get_grid(dh)))
    end

    # Loop over the subdofhandlers
    for sdh in dh.subdofhandlers
        # Check if this sdh contains this field, otherwise continue to the next
        field_idx = Ferrite._find_field(sdh, fieldname)
        field_idx === nothing && continue

        ip = Ferrite.getfieldinterpolation(sdh, field_idx)
        _evaluate_at_grid_nodes_iga!(data, sdh, u, ip, dh.field_offsets[field_idx])
    end
    return data
end

function _evaluate_at_grid_nodes_iga!(
        data::Union{Vector, Matrix}, sdh::SubDofHandler,
        u::AbstractVector{T}, ip::IGAInterpolation, offset::Int
    ) where {T}
    uniformGrid = parameterSpaceGrid(sdh.dh.grid)

    for (nodeid, node) in enumerate(getnodes(uniformGrid))
        val = interpolate(_maybeDeref(ip.basis), u, node.x; offset = offset)
        if data isa Matrix # VTK
            # data[1:length(val), nodeid] .= val
            # data[(length(val) + 1):end, nodeid] .= 0 # purge the NaN
            dataview = @view data[:, nodeid]
            fill!(dataview, 0) # purge the NaN
            Ferrite.toparaview!(dataview, val)
        else
            data[nodeid] = val
        end
    end
    return data
end

function _evaluate_at_grid_nodes_iga!(
        data::Union{Vector, Matrix}, sdh::SubDofHandler,
        u::AbstractVector{T}, ip::VectorizedInterpolation{vdim, shape, order, <:IGAInterpolation},
        offset::Int
    ) where {T, vdim, shape, order}
    uniformGrid = parameterSpaceGrid(sdh.dh.grid)
    basis = _maybeDeref(ip.ip.basis)

    shape_values_gs = gsMatrix()
    actives_gs = gsMatrix{Int32}()
    for (nodeid, node) in enumerate(getnodes(uniformGrid))
        x = Vector(node.x)
        active!(basis, x, actives_gs)
        TinyGismo.eval!(basis, x, shape_values_gs)
        shape_values = toVector(shape_values_gs)
        actives = toVector(actives_gs)

        # Interleaved dof numbering (u1x, u1y, u2x, u2y, ...), see the IGADofHandler `close!`
        val = zero(Vec{vdim, T})
        for k in eachindex(shape_values)
            p = actives[k]
            comps = ntuple(c -> u[offset + vdim * (p - 1) + c], vdim)
            val += shape_values[k] * Vec{vdim, T}(comps)
        end

        if data isa Matrix # VTK
            dataview = @view data[:, nodeid]
            fill!(dataview, 0) # purge the NaN
            Ferrite.toparaview!(dataview, val)
        else
            data[nodeid] = val
        end
    end
    return data
end
