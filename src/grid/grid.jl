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

    kvs = map(FerriteGismo.KnotSpanWrapper{rdim}, knotSpans(basis))

    actives = gsMatrix{Int32}()
    for kv in kvs
        active!(basis, Vector(kv.lower), actives)
        push!(nodes, IGACell{rdim}(toVector(actives)))
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
    return Int(numElements(TinyGismo.basis(grid.geometry), dir))
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

"""
    parameterSpaceGrid(grid::IGAGrid)

Build a standard Ferrite grid of the parameter (knot) space of `grid`, using one
`Line`/`Quadrilateral` cell per knot span. The resulting grid has a regular Lagrange
mesh whose nodes correspond to the knot-span corners in parameter space and is used
internally for VTK export and for evaluating solutions at grid nodes.
"""
function parameterSpaceGrid(grid::IGAGrid{sdim, 2}) where {sdim}
    return generate_grid(
        Quadrilateral, numElements(grid), grid.knotSpans[1].lower, grid.knotSpans[end].upper
    )
end

function parameterSpaceGrid(grid::IGAGrid{sdim, 1}) where {sdim}
    return generate_grid(
        Line, numElements(grid), grid.knotSpans[1].lower, grid.knotSpans[end].upper
    )
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
