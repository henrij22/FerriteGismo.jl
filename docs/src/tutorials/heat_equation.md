# Heat equation

## Introduction

This tutorial solves the stationary heat equation on the unit square using an
Isogeometric Analysis (IGA) discretization. It is the IGA counterpart of the
[Ferrite.jl heat equation tutorial](https://ferrite-fem.github.io/Ferrite.jl/stable/tutorials/heat_equation/),
and you will notice that — apart from constructing the spline geometry, grid and
interpolation — the code is essentially identical to a standard Ferrite program.

The strong form of the problem is

```math
-\nabla \cdot (\nabla u) = 1 \quad \textbf{x} \in \Omega,
```

with homogeneous Dirichlet boundary conditions ``u = 0`` on the whole boundary
``\partial\Omega``. The corresponding weak form reads: find ``u \in \mathbb{U}`` such that

```math
\int_\Omega \nabla \delta u \cdot \nabla u \, \mathrm{d}\Omega =
\int_\Omega \delta u \, \mathrm{d}\Omega \quad \forall \, \delta u \in \mathbb{T},
```

where ``\mathbb{U}`` and ``\mathbb{T}`` are the trial and test function spaces. In IGA these
spaces are spanned by spline (B-spline / NURBS) basis functions instead of the usual
Lagrange polynomials.

## Commented program

We start by loading FerriteGismo, which re-exports both `Ferrite` and `TinyGismo`, together
with `SparseArrays` for the global matrix.

```@example heat
using FerriteGismo, SparseArrays
```

### Geometry and grid

Instead of generating a mesh, we create a spline **geometry** with TinyGismo. Here we start
from a B-spline description of a square of side length `2.0` and refine it uniformly a few
times to obtain a reasonable number of elements. The [`IGAGrid`](@ref) then wraps this
geometry so that it can be used with the Ferrite machinery.

```@example heat
geometry = createBSplineSquare(2.0)
uniformRefine!(geometry, 3)
grid = IGAGrid{2}(geometry)
```

### Interpolation and cell values

The spline basis of the geometry is exposed as a Ferrite interpolation through
[`IGAInterpolation`](@ref). We use it both as the field interpolation and, in a vectorized
form, for the geometry mapping when constructing the `CellValues`. The quadrature rule is
the standard Ferrite `QuadratureRule` for the reference quadrilateral.

```@example heat
ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
qr = QuadratureRule{RefQuadrilateral}(2)
cellvalues = CellValues(qr, ip, ip^2)
```

### Degrees of freedom

Degrees of freedom are distributed with an [`IGADofHandler`](@ref). Just like in Ferrite, we
`add!` a field and `close!` the handler.

```@example heat
dh = IGADofHandler(grid)
add!(dh, :u, ip)
close!(dh)

K = allocate_matrix(dh)
nothing # hide
```

### Boundary conditions

We prescribe homogeneous Dirichlet conditions on all four parametric edges. Because the
standard `ConstraintHandler` interface is not yet fully supported for IGA grids, we use the
[`fixEdge!`](@ref) helper, which adds the boundary control-point dofs directly to the
constraint handler.

```@example heat
ch = ConstraintHandler(dh)
fixEdge!(dh, ch, :left, :u)
fixEdge!(dh, ch, :right, :u)
fixEdge!(dh, ch, :bottom, :u)
fixEdge!(dh, ch, :top, :u)
close!(ch)
```

### Element assembly

The element routine is exactly the same as in the Ferrite tutorial: for each quadrature
point we accumulate the local stiffness matrix `Ke` and load vector `fe`.

```@example heat
function assemble_element!(Ke::Matrix, fe::Vector, cellvalues::CellValues)
    n_basefuncs = getnbasefunctions(cellvalues)
    fill!(Ke, 0)
    fill!(fe, 0)
    # Loop over quadrature points
    for q_point in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q_point)
        # Loop over test shape functions
        for i in 1:n_basefuncs
            δu  = shape_value(cellvalues, q_point, i)
            ∇δu = shape_gradient(cellvalues, q_point, i)
            fe[i] += δu * dΩ
            # Loop over trial shape functions
            for j in 1:n_basefuncs
                ∇u = shape_gradient(cellvalues, q_point, j)
                Ke[i, j] += (∇δu ⋅ ∇u) * dΩ
            end
        end
    end
    return Ke, fe
end
nothing # hide
```

### Global assembly

The global assembly loop uses the Ferrite `CellIterator`, `reinit!` and `assemble!`
functions, which all work transparently on the `IGADofHandler`.

```@example heat
function assemble_global(cellvalues::CellValues, K::SparseMatrixCSC, dh::IGADofHandler)
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)
    f  = zeros(ndofs(dh))
    assembler = start_assemble(K, f)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        assemble_element!(Ke, fe, cellvalues)
        assemble!(assembler, celldofs(cell), Ke, fe)
    end
    return K, f
end

K, f = assemble_global(cellvalues, K, dh)
nothing # hide
```

### Solution of the system

We apply the boundary conditions and solve the linear system for the control-point values
of the solution field.

```@example heat
apply!(K, f, ch)
u = K \ f
maximum(u)
```

### Postprocessing

Finally, the solution is exported to a VTK file. FerriteGismo evaluates the spline field on
the parameter-space grid and writes it out just like Ferrite's `VTKGridFile`.

```@example heat
VTKGridFile("heat_equation", dh) do vtk
    write_solution(vtk, dh, u)
end
nothing # hide
```

The resulting `heat_equation.vtu` file can be opened in e.g. [ParaView](https://www.paraview.org/)
to visualize the solution.
