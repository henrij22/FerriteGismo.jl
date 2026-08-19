# ==============================================================================
# The export mesh
#
# Everything that leaves FerriteGismo for visualization goes through a plain Ferrite
# `Grid` of the knot spans: the spline patch is sampled at a lattice of parametric
# points, and those points become the nodes of ordinary `Line`/`Quadrilateral` cells.
#
# Building that grid explicitly (instead of hooking into Ferrite's internal VTK grid
# construction) is what lets the export use only the public Ferrite API — and it makes
# `subdivision` possible, which samples each knot span several times so that a
# high-order spline field is not drawn as one straight line per element.
# ==============================================================================

"""
    breakpoints(grid::IGAGrid, dir::Integer; subdivision::Int = 1)

Parametric coordinates at which `grid` is sampled along direction `dir`: the knot-span
boundaries of that direction, with `subdivision - 1` equally spaced points inserted inside
every span.

The returned vector is sorted and has `subdivision * numElementsPerDirection(grid, dir) + 1`
entries.
"""
function breakpoints(grid::IGAGrid{sdim, rdim}, dir::Integer; subdivision::Int = 1) where {sdim, rdim}
    @argcheck 1 <= dir <= rdim "Direction $dir out of range for a $rdim-dimensional parameter space"
    @argcheck subdivision >= 1 "subdivision must be at least 1"

    knots = sort!(unique!([ks.lower[dir] for ks in grid.knotSpans]))
    push!(knots, maximum(ks.upper[dir] for ks in grid.knotSpans))

    subdivision == 1 && return knots

    points = Float64[]
    for i in 1:(length(knots) - 1)
        lo, hi = knots[i], knots[i + 1]
        append!(points, range(lo, hi; length = subdivision + 1)[1:(end - 1)])
    end
    push!(points, knots[end])
    return points
end

"""
    exportPoints(grid::IGAGrid; subdivision::Int = 1) -> Vector{Vec{rdim}}

Parametric coordinates of the export mesh nodes of `grid`, ordered lexicographically with
the first parametric direction running fastest.

These are the points at which fields are evaluated for visualization; the physical
positions of the same points are the nodes of [`exportGrid`](@ref).
"""
function exportPoints(grid::IGAGrid{sdim, 1}; subdivision::Int = 1) where {sdim}
    return [Vec{1}((ξ,)) for ξ in breakpoints(grid, 1; subdivision)]
end

function exportPoints(grid::IGAGrid{sdim, 2}; subdivision::Int = 1) where {sdim}
    ξs = breakpoints(grid, 1; subdivision)
    ηs = breakpoints(grid, 2; subdivision)
    return [Vec{2}((ξ, η)) for η in ηs for ξ in ξs]
end

# Cells over the lattice of `exportPoints`, in the same lexicographic ordering
function exportCells(grid::IGAGrid{sdim, 1}; subdivision::Int = 1) where {sdim}
    n = length(breakpoints(grid, 1; subdivision))
    return [Line((i, i + 1)) for i in 1:(n - 1)]
end

function exportCells(grid::IGAGrid{sdim, 2}; subdivision::Int = 1) where {sdim}
    nξ = length(breakpoints(grid, 1; subdivision))
    nη = length(breakpoints(grid, 2; subdivision))
    nodeid(i, j) = (j - 1) * nξ + i
    return [
        Quadrilateral((nodeid(i, j), nodeid(i + 1, j), nodeid(i + 1, j + 1), nodeid(i, j + 1)))
            for j in 1:(nη - 1) for i in 1:(nξ - 1)
    ]
end

"""
    parameterSpaceGrid(grid::IGAGrid; subdivision::Int = 1)

Build a standard Ferrite grid of the parameter (knot) space of `grid`, using one
`Line`/`Quadrilateral` cell per knot span. Its nodes carry the *parametric* coordinates of
the knot-span corners, i.e. the points of [`exportPoints`](@ref).

`subdivision` splits every knot span into that many cells per direction.

See also [`exportGrid`](@ref), which is the same mesh with physical node coordinates.
"""
function parameterSpaceGrid(grid::IGAGrid{sdim, rdim}; subdivision::Int = 1) where {sdim, rdim}
    nodes = [Ferrite.Node(ξ) for ξ in exportPoints(grid; subdivision)]
    return Ferrite.Grid(exportCells(grid; subdivision), nodes)
end

"""
    exportGrid(grid::IGAGrid; subdivision::Int = 1) -> Ferrite.Grid

The mesh that visualization tools are handed in place of the spline patch.

It is an ordinary Ferrite grid of the knot spans of `grid`, with the nodes placed at the
*physical* positions of the knot-span corners.

Since the cells are linear, the picture follows the exact geometry only at the nodes. Use
`subdivision` to sample each knot span several times per direction, which resolves both a
curved geometry and a high-order solution field:

```julia
vtkgrid = exportGrid(grid; subdivision = 3)
```

See also [`parameterSpaceGrid`](@ref) and [`exportPoints`](@ref).
"""
function exportGrid(grid::IGAGrid{sdim, rdim}; subdivision::Int = 1) where {sdim, rdim}
    nodes = [Ferrite.Node(toPhysical(grid, Vector(ξ))) for ξ in exportPoints(grid; subdivision)]
    return Ferrite.Grid(exportCells(grid; subdivision), nodes)
end

"""
    evaluateAtExportNodes(dh::IGADofHandler, u::AbstractVector, fieldname::Symbol; subdivision::Int = 1)

Evaluate the field `fieldname` of the solution vector `u` at the nodes of the export mesh,
i.e. at [`exportPoints`](@ref) of the grid of `dh`.

Returns a `Vector` of scalars or `Vec`s, in the node ordering of [`exportGrid`](@ref), which
is what the VTK writers hand to `Ferrite.write_node_data`.
"""
function evaluateAtExportNodes(
        dh::IGADofHandler, u::AbstractVector, fieldname::Symbol; subdivision::Int = 1
    )
    fieldname ∈ Ferrite.getfieldnames(dh) || error("Field $fieldname not found.")
    sdh = only(dh.subdofhandlers)
    field_index = Ferrite.find_field(sdh, fieldname)
    ip = Ferrite.getfieldinterpolation(dh, (1, field_index))
    offset = dh.field_offsets[field_index]

    return [interpolate(ip, u, Vector(ξ); offset) for ξ in exportPoints(Ferrite.get_grid(dh); subdivision)]
end
