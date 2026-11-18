<h1 align="center">FerriteGismo.jl</h1>

<p align="center">
  <strong>Isogeometric Analysis (IGA) integration for <a href="https://github.com/Ferrite-FEM/Ferrite.jl">Ferrite.jl</a></strong>
</p>

<p align="center">
  <a href="https://henrij22.github.io/FerriteGismo.jl/dev"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Docs dev"></a>
  <a href="https://github.com/henrij22/FerriteGismo.jl/actions/workflows/CI.yml"><img src="https://github.com/henrij22/FerriteGismo.jl/actions/workflows/CI.yml/badge.svg" alt="Build status"></a>
  <a href="https://codecov.io/gh/henrij22/FerriteGismo.jl"><img src="https://codecov.io/gh/henrij22/FerriteGismo.jl/graph/badge.svg?token=IoLnf88JWA" alt="codecov"></a>
  <a href="https://github.com/JuliaTesting/Aqua.jl"><img src="https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg" alt="Aqua QA"></a>
  <a href="https://github.com/fredrikekre/Runic.jl"><img src="https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black" alt="Code style: Runic"></a>
</p>

<p align="center">
  <a href="#features">Features</a> ·
  <a href="#installation">Installation</a> ·
  <a href="https://henrij22.github.io/FerriteGismo.jl/dev">Documentation</a>
</p>

---

FerriteGismo.jl brings Isogeometric Analysis capabilities to the Ferrite finite element framework by interfacing with the  popular [G+Smo](https://github.com/gismo/gismo) (Geometry + Simulation Modules) C++ library.

## ⚠️ Status

This package is under active development and relies on internal Ferrite.jl functionality. Note:

- The `ConstraintHandler` is not currently fully supported and will require future redesign
- API is subject to change as the package matures
- Implementing a custom `IGADofHandler` was also not a good idea and should be worked on in the future

## Features

- **IGA Grids** (`IGAGrid`): IGA grid implementation compatible with `Ferrite.AbstractGrid`
- **IGA Interpolations** (`IGAInterpolation`): B-spline and NURBS basis functions as Ferrite interpolations
- **Geometric Utilities**: Access to knot vectors, degree elevation, knot refinement, and more
- **VTK Export**: Visualization support for IGA results

## Installation

### Prerequisites

This package requires both FerriteGismo.jl and its dependency TinyGismo.jl, which are available through my Julia registry.

### Setup

1. Add my registry:

```julia
using Pkg
pkg"registry add https://github.com/henrij22/JuliaRegistry"
```

1. Install FerriteGismo:

```julia
pkg"add FerriteGismo"
```

## Quick Start

```@example
using FerriteGismo

# Create a B-spline geometry (unit square)
geometry = createBSplineSquare(1.0)

# Optional: refine the geometry
degreeElevate!(geometry, 1)      # Increase polynomial degree
uniformRefine!(geometry, 1)       # h-refinement

# Create an IGA grid
grid = IGAGrid{2}(geometry)

# Create IGA interpolation from the basis
basis = TinyGismo.basis(geometry)
ip = IGAInterpolation{RefQuadrilateral}(basis)

# Setup dof handler (same as in Ferrite)
dh = IGADofHandler(grid)
add!(dh, :u, ip)
close!(dh)

# Create cell values for assembly
qr = QuadratureRule{RefQuadrilateral}(2)
cv = CellValues(qr, ip, ip^2)

# Assembly loop works just like standard Ferrite
for cell in CellIterator(dh)
    reinit!(cv, cell)
    # ... assembly operations
end
```

## Architecture

FerriteGismo.jl interfaces with G+Smo through [TinyGismo.jl](https://github.com/henrij22/TinyGismo.jl), an unofficial Julia wrapper that provides:

- **Native Julia conventions**: 1-based indexing, bang methods (`!`), parametric types
- **Core IGA functionality**: Knot vectors, basis functions (B-spline, NURBS, tensor products)
- **Geometry modeling**: Spline curves, surfaces, and volumes
- **Factory functions**: Standard geometric shapes and primitives

The official Julia bindings [Gismo.jl](https://github.com/gismo/Gismo.jl) are also available. TinyGismo.jl aims to provide a more Julian interface with better documentation and type safety.

## Requirements

- Ferrite.jl 1.4

## Contributing

Contributions are welcome! This package is in active development and there are many opportunities to help:

- Implement missing features (constraint handler, boundary conditions)
- Add examples and tutorials
- Improve documentation
- Performance optimizations
- Bug fixes and testing

## Citation

If you use FerriteGismo.jl in your research, please cite this package:

```bibtex
@software{FerriteGismo,
  author = {Jakob, Henrik},
  title = {FerriteGismo.jl: Isogeometric Analysis for Ferrite.jl},
  url = {https://github.com/henrij22/FerriteGismo.jl},
  version = {0.1.2},
  year = {2026}
}
```

## Related Projects

- [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) - Finite element framework
- [G+Smo](https://github.com/gismo/gismo) - Geometry + Simulation Modules
- [TinyGismo.jl](https://github.com/henrij22/TinyGismo.jl) - Julia wrapper for G+Smo
- [Gismo.jl](https://github.com/gismo/Gismo.jl) - Official Julia bindings for G+Smo

## License

See [LICENSE](LICENSE) file for details.
