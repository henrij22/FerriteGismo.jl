# Postprocessing and export

## Evaluating the solution

The control-point values of a solution vector are not values of the solution — evaluating the
field means summing the basis functions. [`interpolate`](@ref) does that at a **parametric**
point:

```@example post
using FerriteGismo

geometry = createBSplineSquare(1.0)
degreeElevate!(geometry, 1)
uniformRefine!(geometry, 2)
grid = IGAGrid{2}(geometry)
ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))

dh = IGADofHandler(grid)
add!(dh, :u, ip^2)
close!(dh)

u = zeros(ndofs(dh))
u[1:2:end] .= [Ferrite.get_node_coordinate(grid, i)[1] for i in 1:getnnodes(grid)]

interpolate(ip^2, u, [0.5, 0.25])
```

To go the other way — from a physical point to the parametric one — use
`TinyGismo.closestPointTo`. `FerriteGismo.toPhysical(grid, ξ)` maps a parametric point to
physical space.

Ferrite's `evaluate_at_grid_nodes(dh, u, :u)` also works and returns the values at the
corners of the knot spans, which is what the VTK export uses.

## VTK export

```@example post
mktempdir() do dir                                     # hide
VTKGridFile(joinpath(dir, "solution"), dh) do vtk      # hide
    write_solution(vtk, dh, u)
end                                                    # hide
end                                                    # hide
nothing # hide
```

```julia
VTKGridFile("solution", dh) do vtk
    write_solution(vtk, dh, u)
end
```

A spline patch cannot be written to a VTK file, so what is exported is the **export mesh**:
an ordinary Ferrite grid with one linear cell per knot span, whose nodes sit at the exact
physical positions of the knot-span corners. The spline field is evaluated at those nodes.
You can build that mesh yourself — it is a plain `Grid` and useful for any other
visualization route as well:

```@example post
eg = exportGrid(grid)
getncells(eg), getnnodes(eg), typeof(eg)
```

Because the cells are linear, a curved geometry and a high-order solution are drawn as
straight lines between the knot-span corners. `subdivision` samples each span several times
per direction and fixes that; pass the same value when opening the file and when writing
data, since the two have to agree on the nodes:

```julia
VTKGridFile("solution", dh; subdivision = 4) do vtk
    write_solution(vtk, dh, u; subdivision = 4)
end
```

```@example post
getncells(exportGrid(grid; subdivision = 4))
```

The nodes added by subdividing lie on the *exact* geometry, not on a polygonal
approximation of it, so a NURBS arc gets rounder rather than just finer.

## Projecting quadrature-point data

Stresses and other quantities computed at quadrature points are discontinuous data that has
to be projected onto the spline space before it can be written out. `L2Projector` accepts an
`IGAInterpolation` and an `IGAGrid`:

```julia
projector = L2Projector(ip, grid)
σ_nodes = project(projector, σ_qp, qr)     # σ_qp[cellid][q_point]

VTKGridFile("stresses", dh) do vtk
    write_solution(vtk, dh, u)
    write_projection(vtk, projector, σ_nodes, "stress")
end
```

`project` expects the quadrature-point values indexed by **cell id** (not by the index inside
a cell set), in the order in which the quadrature rule visits the points. The projector
builds and factorizes the spline mass matrix once, so reuse it across load steps rather than
constructing a new one each time. Scalars and tensors are both supported — tensor components
are projected independently — and `qr_lhs` defaults to a rule that integrates the mass matrix
of the spline basis exactly.

`evaluate_at_grid_nodes(projector, σ_nodes)` reads the projected field off the export mesh,
and also takes a `subdivision` keyword.

## Reaction forces of a boundary

The residual entries of constrained dofs are the reactions. Summing them over a whole
parametric side gives the resultant force on that side:

```julia
loaded = edgeControlPoints(dh, :left, :u)
Fx = sum(r[2 * (p - 1) + 1] for p in loaded)
```

This is the consistent integral of the traction over the side, and it is the recommended way
to obtain a load-displacement curve from a displacement-controlled analysis. A single
control-point reaction, on the other hand, has no physical meaning on its own.
