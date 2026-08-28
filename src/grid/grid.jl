# Nodes are the indices of active basefunctions
"""
    IGACell{rdim} <: Ferrite.AbstractCell{Ferrite.RefHypercube{rdim}}

Internal cell type of an [`IGAGrid`](@ref), representing one knot span (Bézier element) of
reference dimension `rdim`. Its `nodes` field stores the global indices of the spline basis
functions that are active on the knot span; these play the role of the cell's "nodes" in
the Ferrite sense and are what `celldofs` are derived from.
"""
struct IGACell{rdim} <: Ferrite.AbstractCell{Ferrite.RefHypercube{rdim}}
    nodes::Vector{Int}
end

"""
    IGAGrid{sdim}(geometry::gsGeometry)
    IGAGrid{sdim, rdim}(geometry::gsGeometry)

An Isogeometric Analysis (IGA) grid backed by a G+Smo geometry, implementing the
[`Ferrite.AbstractGrid`](https://ferrite-fem.github.io/Ferrite.jl/stable/reference/grid/)
interface.

Unlike a standard Ferrite grid, the "nodes" of an `IGAGrid` are the spline control
points of `geometry` and each cell corresponds to a non-empty knot span (Bézier element)
of the underlying tensor-product basis. The cells store the indices of the basis functions
that are active on the respective knot span, so that the usual Ferrite assembly machinery
(`CellIterator`, `celldofs`, `reinit!`, …) can be reused.

# Type parameters
- `sdim`: spatial dimension (dimension of the physical space the control points live in).
- `rdim`: reference/parametric dimension of the geometry. Defaults to `sdim` when only a
  single parameter is given.

# Example
```julia
geometry = createBSplineSquare(1.0)
uniformRefine!(geometry, 3)
grid = IGAGrid{2}(geometry)
```

See also [`parameterSpaceGrid`](@ref), [`numElements`](@ref) and
[`numElementsPerDirection`](@ref).
"""
struct IGAGrid{sdim, rdim, G <: gsGeometry} <: Ferrite.AbstractGrid{sdim}
    cells::Vector{IGACell{rdim}}
    nodes::Vector{Ferrite.Node}
    geometry::G
    knotSpans::Vector{KnotSpanWrapper{rdim}}
end

get_sdim(::IGAGrid{sdim}) where {sdim} = sdim::Int
get_rdim(::IGAGrid{sdim, rdim}) where {sdim, rdim} = rdim::Int

IGAGrid{dim}(geometry::gsGeometry) where {dim} = IGAGrid{dim, dim}(geometry)

function IGAGrid{sdim, rdim}(geometry::G) where {G <: gsGeometry, sdim, rdim}
    VecT = typeof(zero(Vec{sdim}))
    cc = Ferrite.Node.(extractCoefs(toMatrix(TinyGismo.coefs(geometry)), VecT))

    nodes = IGACell[] # initialize with undef
    basis = TinyGismo.basis(geometry)

    kvs = _elementSpans(basis, Val(rdim))

    actives = gsMatrix{Int32}()
    for kv in kvs
        push!(nodes, IGACell{rdim}(_activeIn(basis, kv, actives)))
    end

    return IGAGrid{sdim, rdim, G}(nodes, cc, geometry, kvs)
end

# See also spatial_coordinate(cellvalues, qp, cell_coordinates)
function toPhysical(grid::IGAGrid{sdim}, x) where {sdim}
    result = zero(Vec{sdim})

    shape_values_gs = gsMatrix()
    actives_gs = gsMatrix{Int32}()

    active!(TinyGismo.basis(grid.geometry), x, actives_gs)
    TinyGismo.eval!(TinyGismo.basis(grid.geometry), x, shape_values_gs)

    shape_values = toVector(shape_values_gs)
    actives = toVector(actives_gs)

    for (shape_value, idx) in zip(shape_values, actives)
        result += shape_value * grid.nodes[idx].x
    end
    return result
end


"""
    numElementsPerDirection(grid::IGAGrid, dir::Integer)

Return the number of (non-empty) knot spans / Bézier elements of `grid` along the
parametric direction `dir`.
"""
function numElementsPerDirection(grid::IGAGrid{sdim, 2}, dir::Integer) where {sdim}
    @argcheck !isHierarchical(grid) "a hierarchical grid has no per-direction element count: its elements do not form a tensor-product lattice. Use `getncells(grid)` for the total."
    # Counted from the knot spans rather than asked of the basis: for a patch whose knots
    # were inserted non-uniformly, `TinyGismo.numElements(basis, dir)` does not report the
    # number of spans of that direction, and the product over the directions then
    # disagrees with `getncells`.
    return length(breakpoints(grid, dir)) - 1
end

"""
    numElements(grid::IGAGrid)

Return a tuple with the number of elements of `grid` per parametric direction,
e.g. `(nx, ny)` for a two-dimensional grid. The product of the entries equals
`getncells(grid)`.
"""
TinyGismo.numElements(grid::IGAGrid{sdim, 1}) where {sdim} = (getncells(grid),)

function TinyGismo.numElements(grid::IGAGrid{sdim, 2}) where {sdim}
    return (numElementsPerDirection(grid, 1), numElementsPerDirection(grid, 2))
end

Ferrite.geometric_interpolation(::Type{IGACell{2}}) = Lagrange{RefQuadrilateral, 1}()
Ferrite.get_coordinate_type(::IGAGrid{sdim}) where {sdim} = Vec{sdim, Float64}
Ferrite.get_reference_dimension(::IGAGrid{sdim, rdim}) where {sdim, rdim} = rdim

ref_to_param(r, knotSpan::KnotSpanWrapper) = ref_to_param(r, knotSpan.lower, knotSpan.upper)

function ref_to_param(r, lower::AbstractVector, upper::AbstractVector)
    @argcheck length(r) == length(lower) == length(upper)
    mid = (upper .+ lower) ./ 2
    half = (upper .- lower) ./ 2
    return mid .+ half .* r
end

function areaOfKnotSpan(knotSpan::KnotSpanWrapper{2})
    l, u = knotSpan.lower, knotSpan.upper
    return (u[1] - l[1]) * (u[2] - l[2])
end

function areaOfKnotSpan(knotSpan::KnotSpanWrapper{1})
    l, u = knotSpan.lower, knotSpan.upper
    return u[1] - l[1]
end

# Todo physical_to_param

# Reference shape of a `rdim`-dimensional knot span.
_refShape(::Val{1}) = RefLine
_refShape(::Val{2}) = RefQuadrilateral
_refShape(::Val{3}) = RefHexahedron

"""
    hierarchicalSubdomains(grid::IGAGrid) -> Vector{Tuple{Vector{Int}, IGAInterpolation}}

Partition the cells of a hierarchical grid into groups sharing a number of active basis
functions, and give each group its own [`IGAInterpolation`](@ref).

On a hierarchical patch that number varies between elements, while
`Ferrite.getnbasefunctions` is one number per interpolation. Grouping restores a constant
count within each group, which is what a `SubDofHandler` needs. A tensor grid yields a
single group.

Groups come in ascending order of active count, each a `(cellset, ip)` pair for one
`SubDofHandler`:

```julia
dh = IGADofHandler(grid)
for (cells, ip) in hierarchicalSubdomains(grid)
    sdh = SubDofHandler(dh, cells)
    add!(sdh, :u, ip)
end
close!(dh)
```

All groups share one basis, so dofs stay numbered globally over the control points.
"""
function hierarchicalSubdomains(grid::IGAGrid{sdim, rdim}) where {sdim, rdim}
    basis = TinyGismo.basis(grid.geometry)
    shape = _refShape(Val(rdim))

    groups = Dict{Int, Vector{Int}}()
    for (ci, cell) in pairs(grid.cells)
        push!(get!(groups, length(cell.nodes), Int[]), ci)
    end

    return [
        (cells, IGAInterpolation{shape}(basis, n))
            for (n, cells) in sort!(collect(groups); by = first)
    ]
end
