# FerriteGismo.jl

Welcome to the documentation of **FerriteGismo.jl**, a package that brings
**Isogeometric Analysis (IGA)** to the [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl)
finite element framework.

FerriteGismo interfaces with the [G+Smo](https://github.com/gismo/gismo) (*Geometry +
Simulation Modules*) C++ library through the Julia wrapper
[TinyGismo.jl](https://github.com/henrij22/TinyGismo.jl). It provides the building blocks
needed to solve partial differential equations with spline-based (B-spline and NURBS)
discretizations, while reusing as much of the familiar Ferrite workflow as possible:

- an [`IGAGrid`](@ref) that implements the `Ferrite.AbstractGrid` interface on top of a
  G+Smo geometry,
- an [`IGAInterpolation`](@ref) that exposes a spline basis as a standard Ferrite
  interpolation,
- an [`IGADofHandler`](@ref) that distributes degrees of freedom over the spline control
  points,
- VTK export and `L2Projector` support for postprocessing.

Because interpolations, cell values, iterators and assemblers follow the Ferrite
conventions, an assembly loop written for FerriteGismo looks almost identical to a standard
Ferrite one.

!!! warning "Status"
    FerriteGismo is under active development and relies on internal Ferrite.jl
    functionality. The API is subject to change, and the `ConstraintHandler` is only
    partially supported: Dirichlet conditions are given per parametric side, see the
    [Boundary conditions](topic_guides/boundary_conditions.md) topic guide.

## How the documentation is organized

This high-level view of the documentation structure will help you find what you are looking
for:[^1]

- [**Tutorials**](tutorials/index.md) are thoroughly documented, worked examples that guide
  you through solving a problem with FerriteGismo from start to finish.
- [**Topic Guides**](topic_guides/index.md) explain one building block at a time — creating a
  grid, working with the dof handler, boundary conditions, forces and postprocessing.
- [**API Reference**](api_reference.md) contains the technical reference of the
  user-facing functions and types (the documentation strings).
- [**Developer documentation**](developer.md) explains the internals of FerriteGismo and is
  mainly targeted at people who want to contribute to or extend the package.

[^1]: The organization of the documentation loosely follows the
    [Diátaxis Framework](https://diataxis.fr), as does the
    [Ferrite.jl documentation](https://ferrite-fem.github.io/Ferrite.jl/stable/).

## Installation

FerriteGismo and its dependency TinyGismo are available through a custom Julia registry.
Add the registry and then install the package:

```julia
using Pkg
pkg"registry add https://github.com/henrij22/JuliaRegistry"
pkg"add FerriteGismo"
```

## Quick start

```julia
using FerriteGismo

# Create and refine a B-spline geometry (unit square), see TinyGismo.jl
geometry = createBSplineSquare(1.0)
degreeElevate!(geometry, 1)   # increase the polynomial degree
uniformRefine!(geometry, 3)   # h-refinement

# Build the IGA grid, interpolation, and cell values
grid = IGAGrid{2}(geometry)
ip   = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
qr   = QuadratureRule{RefQuadrilateral}(2)
cv   = CellValues(qr, ip, ip^2)

# Distribute the degrees of freedom (just like in Ferrite)
dh = IGADofHandler(grid)
add!(dh, :u, ip)
close!(dh)

# The assembly loop works exactly like a standard Ferrite one
for cell in CellIterator(dh)
    reinit!(cv, cell)
    # ... element assembly
end
```

A complete, worked version of this workflow is available in the
[Heat equation](tutorials/heat_equation.md) tutorial.

## Related packages

- [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) — the finite element framework
  FerriteGismo builds on.
- [G+Smo](https://github.com/gismo/gismo) — the C++ Geometry + Simulation Modules library.
- [TinyGismo.jl](https://github.com/henrij22/TinyGismo.jl) — the Julian wrapper for G+Smo
  used by FerriteGismo.
