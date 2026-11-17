# Boundary conditions

Dirichlet conditions are the one place where IGA cannot simply reuse the Ferrite machinery.
Ferrite derives the constrained dofs from the nodes of a facet, using the fact that a
Lagrange shape function is one at its own node and zero at all others. Spline basis functions
are not interpolatory, and a control point is not a point of the geometry, so FerriteGismo
prescribes conditions on the **control points of a parametric side** instead.

Three functions cover this, all taking a side as `:left`, `:right`, `:bottom` or `:top` (or
the corresponding G+Smo boundary index):

- [`fixEdge!`](@ref) — homogeneous (zero) conditions,
- [`prescribeEdge!`](@ref) — inhomogeneous, possibly time/load-factor dependent conditions,
- [`fixEdgeTangent!`](@ref) — clamped edges, which also fix the slope across the side.

```@example bc
using FerriteGismo

geometry = createBSplineSquare(1.0)
degreeElevate!(geometry, 1)
uniformRefine!(geometry, 2)
grid = IGAGrid{2}(geometry)
ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))

dh = IGADofHandler(grid)
add!(dh, :u, ip^2)
close!(dh)

ch = ConstraintHandler(dh)
fixEdge!(dh, ch, :left, :u)                        # clamp the whole left edge
fixEdge!(dh, ch, :bottom, :u; components = 2)      # symmetry: only u_y
prescribeEdge!(dh, ch, :right, :u, (x, t) -> 0.1t; components = 1)
close!(ch)

length(ch.prescribed_dofs)
```

Afterwards the handler is an ordinary Ferrite `ConstraintHandler`: `update!(ch, t)`
recomputes the prescribed values for a new time or load factor, and `apply!`/`apply_zero!`
work on the matrix, the right-hand side and the solution vector as usual.

```@example bc
Ferrite.update!(ch, 2.0)
u = zeros(ndofs(dh))
apply!(u, ch)
interpolate(ip^2, u, [1.0, 0.5])[1]   # u_x on the right edge
```

## What "prescribed on the control points" means

`prescribeEdge!` evaluates `f` at the *control-point coordinates* and writes the result into
the corresponding dofs. Since the spline basis is a partition of unity, a value that is
**constant along the edge** — the usual displacement-control case — is reproduced exactly
everywhere on that edge, as the evaluation above shows.

A value that varies along the edge is only interpolated at the control points, not
``L^2``-projected onto it, so it is met approximately. 

## Which dofs are constrained

[`edgeControlPoints`](@ref) returns the control-point indices of a side. Together with the
interleaved numbering (`vdim * (p - 1) + c`, see
[Working with the DofHandler](dofhandler.md)) it is the way to address those dofs by hand,
for example to extract the resultant reaction force of a loaded edge from a residual vector
`r`:

```@example bc
loaded = edgeControlPoints(dh, :right, :u)
reactionDofs = [2 * (p - 1) + 1 for p in loaded]

r = zeros(ndofs(dh))   # in practice: the assembled residual
Fx = sum(@view r[reactionDofs])
```

Summing the reactions of a whole side is the consistent integral of the traction over that
side, i.e. the IGA counterpart of summing the nodal reactions of a facet set.

## Clamped edges

For plates, beams and rotation-free shells, clamping also fixes the slope across a side, which
depends on the row of control points next to it. [`fixEdgeTangent!`](@ref) fixes that row to
zero; together with `fixEdge!` this clamps the edge:

```@example bc
ch = ConstraintHandler(dh)
fixEdge!(dh, ch, :left, :u)
fixEdgeTangent!(dh, ch, :left, :u)   # u = 0 and ∂u/∂n = 0 on the left edge
close!(ch)

length(ch.prescribed_dofs)
```

On a shell patch (`IGAGrid{3, 2}`), the tuple form only ties the normal displacement of the row
next to the side to the edge, `n ⋅ (d₁ - d₀) = 0`, and leaves the tangential motion free. These
are affine constraints, so the matrix has to come from `allocate_matrix(dh, ch)`:

```@example bc
coordinates = toMatrix(coefs(geometry))
tilted = [coordinates 0.5 .* coordinates[:, 1]]   # the plane z = x/2
shellGrid = IGAGrid{3, 2}(TensorBSpline{2}(TinyGismo.basis(geometry), tilted))

shellDh = IGADofHandler(shellGrid)
add!(shellDh, :u, ip^3)
close!(shellDh)

shellCh = ConstraintHandler(shellDh)
fixEdge!(shellDh, shellCh, :left, :u)
fixEdgeTangent!(shellDh, shellCh, :left, (:u,))   # or three scalar fields (:ux, :uy, :uz)
close!(shellCh)
K = allocate_matrix(shellDh, shellCh)

count(!isnothing, shellCh.dofcoefficients)
```

Where two clamped sides meet, both edges have to be held and their prescribed motions must
agree along the normal; otherwise only one of the two conditions is enforced. Hierarchical
(THB) patches are not supported.

## Restrictions

- Conditions can only be applied to a whole parametric side, not to a part of it or to a
  single point.
- A `Dirichlet` handed to `add!` must be given on a set of control-point indices; passing a
  facet set is an error, since facets carry no dof information on an IGA grid.
- Affine constraints need the matrix from `allocate_matrix(dh, ch)`.
- `ProjectedDirichlet` is not supported.
