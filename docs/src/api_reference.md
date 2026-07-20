# API Reference

This page documents the user-facing API of FerriteGismo. Most of the surrounding
functionality — assembling, iterating over cells, solving, quadrature rules, VTK export,
etc. — comes directly from [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/stable/),
which is re-exported by FerriteGismo, so please refer to the Ferrite documentation for
those parts.

```@contents
Pages = ["api_reference.md"]
Depth = 2:2
```

## Grid

An [`IGAGrid`](@ref) wraps a G+Smo geometry and implements the `Ferrite.AbstractGrid`
interface. Its cells correspond to knot spans (Bézier elements) and its "nodes" are the
spline control points.

```@docs
IGAGrid
numElements
numElementsPerDirection
parameterSpaceGrid
```

## Interpolations

```@docs
IGAInterpolation
```

## Degrees of freedom

```@docs
IGADofHandler
```

## Boundary conditions

```@docs
fixEdge!
```

## Utilities

```@docs
interpolate
```
