# ==============================================================================
# Exact VTK export as rational Bézier cells
#
# On one knot span every active spline function is a polynomial, hence a combination of the
# Bernstein polynomials of that span, N = Cᵉ B. The extraction operator Cᵉ is computed from
# G+Smo evaluations, which works alike for B-spline, NURBS and hierarchical bases. Geometry
# and fields are written with one common degree per direction, the highest of all bases.
# ==============================================================================

# Bernstein polynomial `i` (0-based) of degree `p` on [0, 1]
_bernstein(p::Int, i::Int, t::Real) = binomial(p, i) * t^i * (1 - t)^(p - i)
_bernstein(degrees::NTuple{rdim, Int}, idx::NTuple{rdim, Int}, t) where {rdim} =
    prod(d -> _bernstein(degrees[d], idx[d], t[d]), 1:rdim)

# Index tuples of the tensor-product Bernstein basis, lexicographic with direction 1 fastest
_bernsteinIndices(degrees::NTuple{1, Int}) = [(i,) for i in 0:degrees[1]]
_bernsteinIndices(degrees::NTuple{2, Int}) = [(i, j) for j in 0:degrees[2] for i in 0:degrees[1]]

#=
Data shared by all knot spans: the Bernstein degrees, one sample point per Bernstein
polynomial and the inverse of the Bernstein values at these points. The points lie strictly
inside the span, since G+Smo evaluates the neighbouring span on its boundary.
=#
struct BezierExtraction{rdim}
    degrees::NTuple{rdim, Int}
    points::Vector{Vec{rdim, Float64}}
    Binv::Matrix{Float64}
end

function BezierExtraction(degrees::NTuple{rdim, Int}) where {rdim}
    indices = _bernsteinIndices(degrees)
    points = [Vec{rdim}(d -> (idx[d] + 1) / (degrees[d] + 2)) for idx in indices]
    B = [_bernstein(degrees, idx, t) for idx in indices, t in points]
    return BezierExtraction{rdim}(degrees, points, inv(B))
end

_degrees(basis, ::Val{1}) = (Int(TinyGismo.degree(basis)),)
_degrees(basis, ::Val{rdim}) where {rdim} = ntuple(d -> Int(TinyGismo.degree(basis, d)), rdim)

# Weights of a rational basis, `nothing` for a polynomial one
_rationalWeights(basis) = _weightsOf(_maybeDeref(basis))
_weightsOf(basis::Union{TinyGismo.NurbsBasis, TinyGismo.TensorNurbsBasis}) = toVector(TinyGismo.weights(basis))
_weightsOf(_) = nothing

#=
The active functions of `basis` on `span` and their extraction operator Cᵉ, one row per active
function. A rational basis is reduced to its polynomial B-splines first:
Nᵢ = Rᵢ / (wᵢ Σⱼ Rⱼ / wⱼ), since Σⱼ Rⱼ / wⱼ = 1 / W.
=#
function _extractionOperator(ext::BezierExtraction{rdim}, basis, weights, span::KnotSpanWrapper{rdim}) where {rdim}
    actives = _activeIn(basis, span)
    N = zeros(length(actives), length(ext.points))
    values = gsMatrix()
    for (k, t) in pairs(ext.points)
        ξ = [span.lower[d] + t[d] * (span.upper[d] - span.lower[d]) for d in 1:rdim]
        eval!(basis, ξ, values)
        R = toVector(values)
        if isnothing(weights)
            N[:, k] .= R
        else
            w = weights[actives]
            N[:, k] .= R ./ (w .* sum(R ./ w))
        end
    end
    return actives, N * ext.Binv
end

_vtkCellType(::Val{1}) = WriteVTK.VTKCellTypes.VTK_BEZIER_CURVE
_vtkCellType(::Val{2}) = WriteVTK.VTKCellTypes.VTK_BEZIER_QUADRILATERAL

# Lexicographic index of each Bernstein node, in VTK's node order: vertices, edges, interior
_vtkOrdering((p,)::NTuple{1, Int}) = [1; p + 1; 2:p]

function _vtkOrdering((p, q)::NTuple{2, Int})
    lex(i, j) = j * (p + 1) + i + 1
    vertices = [lex(0, 0), lex(p, 0), lex(p, q), lex(0, q)]
    # VTK numbers every edge in the direction of increasing parameter
    edges = [
        [lex(i, 0) for i in 1:(p - 1)];
        [lex(p, j) for j in 1:(q - 1)];
        [lex(i, q) for i in 1:(p - 1)];
        [lex(0, j) for j in 1:(q - 1)]
    ]
    interior = [lex(i, j) for j in 1:(q - 1) for i in 1:(p - 1)]
    return [vertices; edges; interior]
end

"""
    VTKBezierFile(filename::String, dh::IGADofHandler; kwargs...)

Open a VTK file holding the patch of `dh` as rational Bézier cells, one per knot span. Unlike
the export mesh of `VTKGridFile`, these represent the geometry and the spline fields exactly;
Paraview subdivides them for display ("Nonlinear Subdivision Level").

Geometry and fields share one degree per direction, the highest of the geometry and all
fields of `dh`. Write data with `write_solution`, `write_projection` and `write_cell_data`,
and add the file to a Paraview collection with `pvd[t] = vtk`. `kwargs` are passed to
`WriteVTK.vtk_grid`.

```julia
VTKBezierFile("solution", dh) do vtk
    write_solution(vtk, dh, u)
end
```

!!! note "Rational fields"
    Paraview interpolates the fields with the rational weights of the geometry. A field is
    therefore written exactly if its basis has the same weight function as the geometry,
    e.g. the geometry basis itself or a refined or degree-elevated version of it, and
    an `ArgumentError` is thrown otherwise, e.g. for a B-spline field on a NURBS geometry.
"""
struct VTKBezierFile{rdim, VTK <: WriteVTK.DatasetFile, G <: IGAGrid}
    vtk::VTK
    grid::G
    extraction::BezierExtraction{rdim}
    weights::Vector{Vector{Float64}} # Bézier weights of the geometry per cell
end

function VTKBezierFile(filename::String, dh::IGADofHandler; kwargs...)
    grid = Ferrite.get_grid(dh)
    rdim = get_rdim(grid)
    degrees = _degrees(TinyGismo.basis(grid.geometry), Val(rdim))
    for name in Ferrite.getfieldnames(dh)
        basis = Ferrite.get_base_interpolation(_fieldInterpolation(dh, name)).basis
        degrees = max.(degrees, _degrees(basis, Val(rdim)))
    end
    return _bezierFile(filename, grid, degrees; kwargs...)
end

function VTKBezierFile(f::Function, args...; kwargs...)
    vtk = VTKBezierFile(args...; kwargs...)
    try
        f(vtk)
    finally
        close(vtk)
    end
    return vtk
end

# Bézier points and weights of the geometry on `span`, in lexicographic order
function _bezierGeometry(ext::BezierExtraction, grid::IGAGrid, basis, weights, span::KnotSpanWrapper)
    actives, C = _extractionOperator(ext, basis, weights, span)
    w = isnothing(weights) ? ones(length(actives)) : weights[actives]
    wb = C' * w
    xb = (C' * [w[k] * grid.nodes[a].x for (k, a) in pairs(actives)]) ./ wb
    return xb, wb
end

function _bezierFile(filename::String, grid::IGAGrid{sdim, rdim}, degrees::NTuple{rdim, Int}; kwargs...) where {sdim, rdim}
    ext = BezierExtraction(degrees)
    basis = TinyGismo.basis(grid.geometry)
    weights = _rationalWeights(basis)

    nb = length(ext.points)
    ordering = _vtkOrdering(degrees)
    points = zeros(3, nb * getncells(grid))
    cellWeights = Vector{Vector{Float64}}(undef, getncells(grid))
    cells = WriteVTK.MeshCell[]

    for (e, span) in ((e, c.span) for (e, c) in pairs(grid.cells))
        xb, wb = _bezierGeometry(ext, grid, basis, weights, span)
        offset = nb * (e - 1)
        for k in 1:nb
            points[1:sdim, offset + k] .= xb[k]
        end
        cellWeights[e] = wb
        push!(cells, WriteVTK.MeshCell(_vtkCellType(Val(rdim)), offset .+ ordering))
    end

    vtk = WriteVTK.vtk_grid(filename, points, cells; kwargs...)
    vtk["RationalWeights", WriteVTK.VTKPointData()] = reduce(vcat, cellWeights)
    vtk[WriteVTK.VTKPointData()] = "RationalWeights" => "RationalWeights"
    vtk["HigherOrderDegrees", WriteVTK.VTKCellData()] = repeat(Float64[degrees...; zeros(3 - rdim)], 1, getncells(grid))
    vtk[WriteVTK.VTKCellData()] = "HigherOrderDegrees" => "HigherOrderDegrees"
    return VTKBezierFile(vtk, grid, ext, cellWeights)
end

Base.close(vtk::VTKBezierFile) = (WriteVTK.vtk_save(vtk.vtk); vtk)

function Base.setindex!(pvd::WriteVTK.CollectionFile, vtk::VTKBezierFile, t::Real)
    return WriteVTK.collection_add_timestep(pvd, vtk.vtk, t)
end

# Bézier coefficients of the spline field with coefficients `coef(a)`, `a` the global index of
# a basis function, in the point order of the file
_bezierField(vtk::VTKBezierFile, basis, coef) = _bezierField(vtk.extraction, vtk.grid, vtk.weights, basis, coef)

function _bezierField(ext::BezierExtraction{rdim}, grid::IGAGrid, cellWeights, basis, coef::F) where {rdim, F}
    all(_degrees(basis, Val(rdim)) .<= ext.degrees) ||
        throw(ArgumentError("the field has a higher degree than the file, $(ext.degrees)"))
    weights = _rationalWeights(basis)

    data = typeof(coef(1))[]
    for (e, span) in ((e, c.span) for (e, c) in pairs(grid.cells))
        actives, C = _extractionOperator(ext, basis, weights, span)
        w = isnothing(weights) ? ones(length(actives)) : weights[actives]
        all(isapprox.(C' * w, cellWeights[e]; rtol = 1.0e-8)) || throw(
            ArgumentError(
                "the field has another rational weight function than the geometry, so it cannot be " *
                    "written exactly. Use the geometry basis (or a refinement of it) or `VTKGridFile`."
            )
        )
        append!(data, (C' * [w[k] * coef(a) for (k, a) in pairs(actives)]) ./ cellWeights[e])
    end
    return data
end

"""
    write_solution(vtk::VTKBezierFile, dh::IGADofHandler, u::AbstractVector, suffix = "")

Write every field of the solution vector `u` to `vtk`, as exact Bézier coefficients.
"""
function Ferrite.write_solution(vtk::VTKBezierFile, dh::IGADofHandler, u::AbstractVector, suffix = "")
    for name in Ferrite.getfieldnames(dh)
        ip = _fieldInterpolation(dh, name)
        offset = fieldOffset(dh, name)
        n = Ferrite.n_components(ip)
        coef(a) = n == 1 ? u[offset + a] : Vec{n}(c -> u[offset + n * (a - 1) + c])
        data = _bezierField(vtk, Ferrite.get_base_interpolation(ip).basis, coef)
        Ferrite._vtk_write_node_data(vtk.vtk, data, string(name, suffix))
    end
    return vtk
end

"""
    write_projection(vtk::VTKBezierFile, proj::L2ProjectorIGA, vals, name)

Write the projected values `vals` (the output of `project`) to `vtk` under `name`.
"""
function Ferrite.write_projection(vtk::VTKBezierFile, proj::L2ProjectorIGA, vals, name)
    data = _bezierField(vtk, proj.ip.basis, a -> vals[a])
    Ferrite._vtk_write_node_data(vtk.vtk, data, name)
    return vtk
end

function Ferrite.write_cell_data(vtk::VTKBezierFile, celldata, name)
    WriteVTK.vtk_cell_data(vtk.vtk, celldata, name)
    return vtk
end
