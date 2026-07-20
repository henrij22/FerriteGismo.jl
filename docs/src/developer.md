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

### Parameter-space grid

For visualization and for evaluating solutions at grid nodes, FerriteGismo builds a regular
Ferrite grid over the knot spans in parameter space via
[`parameterSpaceGrid`](@ref). Each knot span becomes one `Line`/`Quadrilateral` cell, and
its nodes are mapped to physical space with `FerriteGismo.toPhysical` when writing VTK
output.

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

Because the standard `ConstraintHandler` / `Dirichlet` workflow relies on facet sets that do
not exist for `IGAGrid`s, homogeneous Dirichlet conditions along whole parametric edges are
applied with the [`fixEdge!`](@ref) helper. It queries the boundary control points of a
field from G+Smo and calls `Ferrite.add_prescribed_dof!` for each of them. Extending this to
inhomogeneous and more general boundary conditions is future work.

## Postprocessing and L2 projection

VTK export reuses the [parameter-space grid](@ref "Parameter-space grid"): the solution is
evaluated at its nodes and written with Ferrite's VTK backend. For projecting
quadrature-point data (e.g. stresses) onto a continuous spline field, FerriteGismo provides
an IGA-specific projector.

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
