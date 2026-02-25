# Nodes are the indices of active basefunctions
struct IGACell{rdim} <: Ferrite.AbstractCell{Ferrite.RefHypercube{rdim}}
    nodes::Vector{Int}
end

struct IGAGrid{sdim, rdim, G <: gsGeometry} <: Ferrite.AbstractGrid{sdim}
    cells::Vector{IGACell{rdim}}
    nodes::Vector{Ferrite.Node}
    geometry::G
    knotSpans::Vector{KnotSpanWrapper{rdim}}
end

IGAGrid{dim}(geometry::gsGeometry) where {dim} = IGAGrid{dim, dim}(geometry)

function IGAGrid{sdim, rdim}(geometry::G) where {G <: gsGeometry, sdim, rdim}
    VecT = typeof(zero(Vec{sdim}))
    cc = Ferrite.Node.(extractCoefs(toMatrix(TinyGismo.coefs(geometry)), VecT))

    nodes = IGACell[] # initialize with undef
    basis = TinyGismo.basis(geometry)

    kvs = map(FerriteGismo.KnotSpanWrapper{rdim}, knotSpans(basis))

    # kvs = KnotSpan{rdim}[]
    actives = gsMatrix{Int32}()
    for kv in kvs
        # activeInCell = copyMatrix(actives(Gismo.basis(geometry), eles[:, 2i - 1]))[:, 1] .+ 1
        active!(basis, Vector(kv.lower), actives)
        push!(nodes, IGACell{rdim}(toVector(actives)))
        # push!(knotSpans, KnotSpan{rdim}(Vec{rdim}(eles[:, 2i - 1]), Vec{rdim}(eles[:, 2i])))
    end

    return IGAGrid{sdim, rdim, G}(nodes, cc, geometry, kvs)
end

# See also spatial_coordinate(cellvalues, qp, cell_coordinates)
function toPhysical(grid::IGAGrid{sdim}, x) where {sdim}
    result = zero(Vec{sdim})

    shape_values = gsMatrix()
    actives = gsMatrix{Int32}()

    active!(TinyGismo.basis(grid.geometry), x, actives)
    TinyGismo.eval!(TinyGismo.basis(grid.geometry), x, shape_values)

    shape_values = toVector(shape_values)
    actives = toVector(actives)

    for i in eachindex(shape_values)
        result += shape_values[i] * grid.nodes[actives[i]].x
    end
    return result
end


function numElementsPerDirection(grid::IGAGrid{sdim, 2}, dir::Integer) where {sdim}
    return Int(numElements(TinyGismo.basis(grid.geometry), dir))
end

TinyGismo.numElements(grid::IGAGrid{sdim, 1}) where {sdim} = (getncells(grid),)

function TinyGismo.numElements(grid::IGAGrid{sdim, 2}) where {sdim}
    return (numElementsPerDirection(grid, 1), numElementsPerDirection(grid, 2))
end

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
