mutable struct IGADofHandler{dim, G <: IGAGrid{dim}} <: Ferrite.AbstractDofHandler
    dh::DofHandler{dim, G}
    field_offsets::Vector{Int}
end

# Access to wrapped dofhandler components
function Base.getproperty(dh::IGADofHandler, sym::Symbol)
    return (sym === :dh || sym == :field_offsets) ? getfield(dh, sym) : getproperty(dh.dh, sym)
end
function Base.setproperty!(dh::IGADofHandler, sym::Symbol, v)
    return setproperty!(dh.dh, sym, v)
end

function IGADofHandler(grid::G) where {dim, G <: IGAGrid{dim}}
    return IGADofHandler{dim, G}(DofHandler(grid), Int[])
end

Ferrite.ndofs_per_cell(dh::IGADofHandler) = ndofs_per_cell(dh.dh)
Ferrite.ndofs_per_cell(dh::IGADofHandler, ::Int) = ndofs_per_cell(dh.dh)

function Ferrite.add!(dh::IGADofHandler, name::Symbol, ip::Interpolation)
    return add!(dh.dh, name, ip)
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
                activeInCell = toVector(actives) .+ current_dof_counter

                # Interweave dofs u1x, u1y, u1z, u2x, ...
                for j in 1:length(activeInCell)
                    for c in 1:ncomp
                        push!(v[i], ncomp * activeInCell[j] + c - ncomp)
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

Ferrite.allocate_matrix(dh::IGADofHandler) = Ferrite.allocate_matrix(dh.dh)
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
    field_idx = Ferrite.find_field(dh.dh, fieldname)
    ip = Ferrite.getfieldinterpolation(dh.dh, field_idx)
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
        val = interpolate(ip.basis, u, node.x; offset = offset)
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
