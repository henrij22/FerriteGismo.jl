# Linear elasticity

## Introduction

This tutorial solves small-strain linear elasticity with an Isogeometric Analysis (IGA)
discretization. It is the IGA counterpart of the
[Ferrite.jl linear elasticity tutorial](https://ferrite-fem.github.io/Ferrite.jl/stable/tutorials/linear_elasticity/),
and we refer to that tutorial for the theory (weak form, the elastic stiffness tensor, and
the meaning of the assembly expressions). As in the [Heat equation](heat_equation.md)
tutorial, the code is almost identical to a standard Ferrite program — only the geometry,
grid and interpolation are IGA-specific.

We consider a unit square that is clamped along its left edge (a cantilever) and loaded by a
constant body force pointing downwards. Compared to the Ferrite tutorial, the Neumann
traction is replaced by a body force, because facet-based integration
(`FacetValues` / `FacetIterator`) is not yet available for `IGAGrid`s.

## Commented program

```@example elasticity
using FerriteGismo, SparseArrays
```

### Geometry, grid and cell values

We start from a B-spline square, elevate its degree to obtain a quadratic basis, and refine
it uniformly. The displacement is a **vector** field, so the field interpolation is the
vectorized [`IGAInterpolation`](@ref) `ip^2`.

```@example elasticity
geometry = createBSplineSquare(1.0)
degreeElevate!(geometry, 1)   # quadratic basis
uniformRefine!(geometry, 3)
grid = IGAGrid{2}(geometry)

ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
qr = QuadratureRule{RefQuadrilateral}(3)
cellvalues = CellValues(qr, ip^2, ip^2)

dh = IGADofHandler(grid)
add!(dh, :u, ip^2)
close!(dh)
nothing # hide
```

### Boundary conditions

The whole left edge is clamped, i.e. both displacement components are fixed there. This is
done with [`fixEdge!`](@ref); for a vector field it constrains all components by default
(use the `components` keyword to fix only some of them).

```@example elasticity
ch = ConstraintHandler(dh)
fixEdge!(dh, ch, :left, :u)
close!(ch)
```

### Material

The linear-elastic stiffness tensor is assembled from the shear and bulk moduli exactly as
in the Ferrite tutorial.

```@example elasticity
Emod = 200.0e3  # Young's modulus
ν    = 0.3      # Poisson's ratio

Gmod = Emod / (2 * (1 + ν))    # shear modulus
Kmod = Emod / (3 * (1 - 2ν))   # bulk modulus

C = gradient(ϵ -> 2 * Gmod * dev(ϵ) + 3 * Kmod * vol(ϵ), zero(SymmetricTensor{2, 2}))
nothing # hide
```

### Element and global assembly

The element routine adds the elastic stiffness `∇Nᵢ ⊡ C ⊡ ∇ˢʸᵐNⱼ` and the body-force
contribution `b ⋅ Nᵢ`.

```@example elasticity
b = Vec{2}((0.0, -1.0e3))   # constant downward body force

function assemble_cell!(ke, fe, cellvalues, C, b)
    for q_point in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:getnbasefunctions(cellvalues)
            Nᵢ  = shape_value(cellvalues, q_point, i)
            ∇Nᵢ = shape_gradient(cellvalues, q_point, i)
            fe[i] += (b ⋅ Nᵢ) * dΩ
            for j in 1:getnbasefunctions(cellvalues)
                ∇ˢʸᵐNⱼ = shape_symmetric_gradient(cellvalues, q_point, j)
                ke[i, j] += (∇Nᵢ ⊡ C ⊡ ∇ˢʸᵐNⱼ) * dΩ
            end
        end
    end
    return ke, fe
end

function assemble_global(cellvalues, K, dh, C, b)
    n_basefuncs = getnbasefunctions(cellvalues)
    ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)
    f  = zeros(ndofs(dh))
    assembler = start_assemble(K, f)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        fill!(ke, 0.0)
        fill!(fe, 0.0)
        assemble_cell!(ke, fe, cellvalues, C, b)
        assemble!(assembler, celldofs(cell), ke, fe)
    end
    return K, f
end

K = allocate_matrix(dh)
K, f = assemble_global(cellvalues, K, dh, C, b)
nothing # hide
```

### Solution and postprocessing

We apply the boundary conditions, solve, and report the maximum (downward) tip deflection.

```@example elasticity
apply!(K, f, ch)
u = K \ f

uy = u[2:2:end]  # vertical displacement components
minimum(uy)
```

Finally the displacement field is exported to VTK for visualization.

```@example elasticity
VTKGridFile("linear_elasticity", dh) do vtk
    write_solution(vtk, dh, u)
end
nothing # hide
```

The resulting `linear_elasticity.vtu` file can be opened in
[ParaView](https://www.paraview.org/); warping the mesh by the displacement field shows the
bending of the clamped square under its own load.
