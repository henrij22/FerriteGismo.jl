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

The exported grid is the *knot span* grid: one linear VTK cell per element, with its corners
placed at the exact physical positions of the knot-span corners. The spline field is
evaluated at those corners, so a high-order or strongly curved solution looks piecewise
linear in ParaView. Refine the geometry (or export a refined evaluation) if you need a
smoother picture — the underlying solution is as smooth as the basis, the export is not.

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
constructing a new one each time.

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
