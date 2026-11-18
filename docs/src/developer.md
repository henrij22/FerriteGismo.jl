# Developer documentation

How FerriteGismo plugs a G+Smo spline discretization into Ferrite, for contributors.

!!! warning
    FerriteGismo depends on a number of **internal** Ferrite functions, which makes it
    sensitive to Ferrite releases. The custom [`IGADofHandler`](@ref) in particular is a
    candidate for a redesign.

## Design overview

Ferrite's assembly loop, iterators, cell values and solvers are reused unchanged; only the
pieces specific to the geometry are replaced:

| Concept              | Ferrite                    | FerriteGismo                          |
|:---------------------|:---------------------------|:--------------------------------------|
| Geometry             | Nodes + Lagrange cells     | G+Smo geometry (control points)       |
| Grid                 | `Grid`                     | [`IGAGrid`](@ref)                     |
| Cell                 | `Quadrilateral`, …         | `IGACell` (one per knot span)         |
| Interpolation        | `Lagrange`, …              | [`IGAInterpolation`](@ref)            |
| Dof distribution     | `DofHandler`               | [`IGADofHandler`](@ref)               |

The backend (knot vectors, basis evaluation, geometry, refinement, …) is provided by
[TinyGismo.jl](https://github.com/henrij22/TinyGismo.jl), a Julian wrapper around the G+Smo
C++ library. FerriteGismo re-exports both `Ferrite` and `TinyGismo`.

## Grid and cells

An [`IGAGrid`](@ref) stores the G+Smo geometry, the control points as Ferrite `Node`s, one
[`FerriteGismo.IGACell`](@ref) per non-empty knot span. Each cell holds its
[`FerriteGismo.KnotSpanWrapper`](@ref) and the indices of the basis functions that are
*active* on that knot span; these take the role of the cell "nodes" and are the basis for
the connectivity used during assembly.

The reference dimension `rdim` (parametric dimension) and spatial dimension `sdim` are kept
as separate type parameters so that embedded geometries (`rdim < sdim`) can be represented,
even though not all code paths support them yet.

```@docs
FerriteGismo.IGACell
FerriteGismo.KnotSpanWrapper
```

### The export mesh

For visualization and for evaluating solutions at grid nodes, FerriteGismo builds an
ordinary Ferrite grid over the knot spans: [`exportGrid`](@ref) places its nodes at the
physical positions of the knot-span corners, [`parameterSpaceGrid`](@ref) is the same mesh
with parametric coordinates, and [`exportPoints`](@ref) are those coordinates as a plain
vector. All three are derived from [`breakpoints`](@ref), which reads the knot-span
boundaries per direction, so the mesh follows the actual (possibly non-uniform) knots.

A `subdivision` keyword samples every knot span several times per direction, which is how a
curved geometry or a high-order field is resolved in the picture.

## Interpolations

[`IGAInterpolation`](@ref) is a `Ferrite.ScalarInterpolation` wrapping a G+Smo basis. Unlike
a Lagrange interpolation, its active basis functions depend on the knot span. The
interpolation stays stateless: `reinit!` maps the quadrature points into the current knot
span (`FerriteGismo.ref_to_param`) before the basis is evaluated through:

- `Ferrite.reference_shape_values!`
- `Ferrite.reference_shape_gradients_and_values!`
- `Ferrite.reference_shape_hessians_gradients_and_values!`

### Geometry mapping

The mapping is Ferrite's `IdentityMapping`, so Jacobian, Hessian (also for embedded
elements) and the push-forward come from Ferrite. Only
`Ferrite.reinit!(::CellValues, ::CellCache{<:Any, <:IGAGrid})` is specialized: it maps the
quadrature points into the knot span, re-evaluates the basis, calls Ferrite's generic
`reinit!` and scales `detJdV` by the knot-span area relative to the reference element.

### Facet values

`Ferrite.reinit!(::FacetValues, ::CellCache, facet_nr)` works like the cell version, with
the weights scaled by the knot-span length along the facet
(`FerriteGismo.facetScaleOfKnotSpan`). Facet sets are computed on demand from the knot-span
corners by `Ferrite.getfacetset(::IGAGrid, name)`.

## Degree-of-freedom handling

The [`IGADofHandler`](@ref) wraps a `Ferrite.DofHandler` and forwards most calls to it, but
owns its `SubDofHandler`s so that `sdh.dh` is the IGA handler. `close!` numbers the dofs by
control point rather than by Lagrange node: each field gets one global block, shared by all
subdofhandlers, with the components interleaved per control point. Several subdofhandlers
are only needed on hierarchical patches, where they group elements by active-function count.

`allocate_matrix` builds the sparsity pattern from the cell dofs itself, since Ferrite's
version is typed on the concrete `DofHandler`. The `CellCache`/`CellIterator` methods in
`iterators.jl` make Ferrite's iteration use the IGA connectivity.

!!! note
    The implementation reaches into a number of `DofHandler` internals. A cleaner approach
    would be to generalize Ferrite's dof distribution rather than reimplement it here.

## Boundary conditions

Neumann conditions need nothing special. Ferrite's `Dirichlet` on a facet set does not work,
since it evaluates the condition at facet dof positions, which spline bases do not have.
Instead conditions are put on the control points of a parametric side:

- [`fixEdge!`](@ref) calls `Ferrite.add_prescribed_dof!` for the boundary control points.
- [`prescribeEdge!`](@ref) and `add!(ch, ::Dirichlet)` on a set of control points store the
  condition like a Ferrite node-set condition, so `update!` computes the values.
- The shell form of [`fixEdgeTangent!`](@ref) collects its conditions until `close!(ch)`,
  where those meeting at a corner are combined into affine constraints.

## Postprocessing and L2 projection

The export path stays on Ferrite's public API: the writers get the
[export mesh](@ref "The export mesh"), a plain Ferrite `Grid`, and the spline fields are
evaluated at its nodes with [`evaluateAtExportNodes`](@ref). FerriteGismo only adds methods of
public functions (`VTKGridFile`, `write_solution`, `write_projection`,
`evaluate_at_grid_nodes`, `L2Projector`, `project`) for its own types.

The exact export, [`VTKBezierFile`](@ref), cannot go through an export mesh: it writes the
Bézier cells with WriteVTK directly. The extraction operator of each knot span is computed
from G+Smo evaluations at points inside the span, which covers B-spline, NURBS and
hierarchical bases alike. For the point data it reuses `Ferrite._vtk_write_node_data`, which
converts tensors to their Paraview components.

The L2 projector assembles and factorizes the spline mass matrix once and solves for each
tensor component. It subtypes the internal `Ferrite.AbstractProjector` so that code
constrained on that type accepts it.

```@docs
FerriteGismo.L2ProjectorIGA
```

## Contributing

Contributions are very welcome. Good starting points are:

- reducing the reliance on internal Ferrite functionality,
- redesigning the [`IGADofHandler`](@ref),
- improving boundary-condition support (a full `ConstraintHandler` integration),
- adding tutorials and examples,
- performance work and additional tests.
