# Developer documentation

This page describes the internals of FerriteGismo. It is aimed at contributors and at users
who want to understand how the package plugs a spline (G+Smo) discretization into the
Ferrite finite element machinery. It is *not* required reading for using the package — see
the [Tutorials](tutorials/index.md) and [API Reference](api_reference.md) for that.

!!! warning
    FerriteGismo currently depends on a number of **internal** Ferrite functions (methods
    that are not part of Ferrite's public API). This makes the package sensitive to Ferrite
    releases and is one of the main areas that will need to be revisited as the package
    matures. Several design decisions — notably the custom [`IGADofHandler`](@ref) — are
    known to be suboptimal and are candidates for a future redesign.

## Design overview

The central idea of FerriteGismo is to reuse Ferrite's assembly loop, iterators, cell
values and solvers unchanged, and to only replace the pieces that are specific to the
geometry description:

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
[`FerriteGismo.IGACell`](@ref) per non-empty knot span, and a cached
[`FerriteGismo.KnotSpanWrapper`](@ref) for each knot span. Each cell holds the indices of
the basis functions that are *active* on its knot span; these take the role of the cell
"nodes" and are the basis for the connectivity used during assembly.

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

[`IGAInterpolation`](@ref) is a `Ferrite.ScalarInterpolation` wrapping a G+Smo basis. The
key difference from a Lagrange interpolation is that the set of active basis functions
depends on the current knot span. To integrate with Ferrite, the interpolation is *mutable*
and carries a `currentElement` field of type
[`FerriteGismo.KnotSpanWrapper`](@ref), which is updated during `reinit!`.

FerriteGismo provides IGA-specific methods for the low-level Ferrite entry points that
evaluate shape functions and their derivatives on the reference element:

- `Ferrite.reference_shape_values!`
- `Ferrite.reference_shape_gradients_and_values!`
- `Ferrite.reference_shape_hessians_gradients_and_values!`

These map the reference quadrature points into the current knot span (via
`FerriteGismo.ref_to_param`) and delegate the actual evaluation to the G+Smo basis.

### Geometry mapping

The mapping from the reference element to physical space is implemented through a dedicated
`FerriteGismo.IGAMapping` type. FerriteGismo specializes `Ferrite.calculate_mapping` and
`Ferrite.apply_mapping!` for this mapping type to compute the Jacobian (and, where
available, the Hessian) and to push the reference gradients forward to physical space. The
determinant of the Jacobian is additionally scaled by the ratio between the knot-span area
and the reference-element area, so that integration is performed correctly over the
parametric domain.

### Facet values

Boundary integration follows the same pattern: FerriteGismo specializes
`Ferrite.reinit!(::FacetValues, ::CellCache, facet_nr)` for IGA interpolations. Because the
facet quadrature points live on the reference cell, the spline basis is re-evaluated for
the active knot span on every `reinit!`, exactly as for `CellValues`. The weighted facet
normal is computed from the parametric Jacobian with `Ferrite.weighted_normal`, and the
quadrature weights are scaled by the knot-span extent along the facet's tangential
direction (`FerriteGismo.facetScaleOfKnotSpan`). Boundary facet sets are computed on demand
by `Ferrite.getfacetset(::IGAGrid, name)` from the knot-span corners, since an `IGAGrid` is
a structured tensor-product grid.

## Degree-of-freedom handling

The [`IGADofHandler`](@ref) wraps a regular `Ferrite.DofHandler` and forwards most property
accesses and interface methods to it. The non-trivial part is `close!`, which distributes
the global dofs according to the globally numbered spline control points rather than
Lagrange nodes. This is handled by the internal `_close_subdofhandler_iga!` routine, which:

1. iterates over the field interpolations of each subdofhandler,
2. queries the active basis functions per knot span from G+Smo,
3. interleaves the components for vector-valued (`VectorizedInterpolation`) fields, and
4. records per-field offsets in `field_offsets` so that fields can later be addressed inside
   the global solution vector.

The custom `CellCache`/`CellIterator`/`celldofs!` methods in `iterators.jl` make sure the
Ferrite iteration protocol picks up the IGA connectivity.

!!! note
    Only a single subdofhandler / cell type is currently supported, and the implementation
    reaches into a number of `DofHandler` internals. A cleaner approach would be to
    generalize Ferrite's dof distribution rather than to reimplement it here.

## Boundary conditions

Neumann boundary conditions use the standard Ferrite `FacetValues` / `FacetIterator`
workflow together with the on-demand facet sets described above.

The standard `ConstraintHandler` / `Dirichlet` workflow is however not yet supported, since
it evaluates the prescribed function at facet dof positions, which spline bases do not
provide. Homogeneous Dirichlet conditions along whole parametric edges are therefore
applied with the [`fixEdge!`](@ref) helper. It queries the boundary control points of a
field from G+Smo and calls `Ferrite.add_prescribed_dof!` for each of them. Extending this to
inhomogeneous and more general boundary conditions is future work.

## Postprocessing and L2 projection

Unlike the rest of the package, the export path deliberately does **not** hook into
Ferrite's internals. A spline patch cannot be written to a VTK file, so instead of
specializing Ferrite's grid-writing machinery, FerriteGismo hands the writers the
[export mesh](@ref "The export mesh") — a plain Ferrite `Grid` they already understand —
and evaluates the spline fields at its nodes with
[`evaluateAtExportNodes`](@ref). Everything is then written through the public
`VTKGridFile`/`VTKHDFGridFile` constructors and `Ferrite.write_node_data`.

What FerriteGismo adds are methods of *public* Ferrite functions dispatched on its own
types — `VTKGridFile`, `VTKHDFGridFile`, `write_solution`, `write_projection`,
`evaluate_at_grid_nodes`, `L2Projector`, `project` — so downstream code that writes
`write_solution(vtk, dh, u)` keeps working unchanged.

The projector for quadrature-point data (e.g. stresses) is likewise self-contained: it
assembles the spline mass matrix with `CellValues`, `allocate_matrix` and
`start_assemble`/`assemble!`, factorizes it once, and solves for the control-point values
of each tensor component. The only internal it still touches is the abstract type
`Ferrite.AbstractProjector`, which it subtypes so that code constrained on that type (as
Ferrite's own projector is) accepts it.

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
