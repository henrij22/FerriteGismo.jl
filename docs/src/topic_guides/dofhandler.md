# Working with the DofHandler

Degrees of freedom in IGA belong to control points, not to element nodes, and a control point
is shared by every knot span in the support of its basis function. The
[`IGADofHandler`](@ref) distributes the dofs accordingly while keeping the Ferrite interface:
fields are added with `add!`, the handler is finalized with `close!`, and `CellIterator`,
`celldofs`, `ndofs` and `allocate_matrix` behave as usual.

```@example dofs
using FerriteGismo

geometry = createBSplineSquare(1.0)
degreeElevate!(geometry, 1)   # quadratic
uniformRefine!(geometry, 2)   # 4 x 4 elements
grid = IGAGrid{2}(geometry)

ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))

dh = IGADofHandler(grid)
add!(dh, :u, ip^2)            # a two-component displacement field
close!(dh)

ndofs(dh), ndofs_per_cell(dh)
```

`ndofs` counts two dofs per control point, while `ndofs_per_cell` counts only the basis
functions active on one knot span — for a quadratic basis that is ``(p+1)^2 = 9`` functions,
i.e. 18 dofs per cell. Neighbouring cells share most of them, which is exactly the
inter-element continuity of the spline space.

## Dof numbering

Within a field the dofs are interleaved per control point:

```
u1x, u1y, u2x, u2y, …
```

so the dofs of control point `p` of a `vdim`-component field are
`offset + vdim * (p - 1) + c` for `c in 1:vdim`, where `offset` is the offset of the field
inside the global vector. That formula is the one thing worth remembering about the
`IGADofHandler`; it is how the boundary conditions, the reaction forces of a boundary and any
hand-written postprocessing address their dofs:

```@example dofs
controlPoint = 7
uy = 2 * (controlPoint - 1) + 2   # global dof of u_y at that control point
```

For several fields the offsets are stored in `dh.field_offsets`, indexed by the field index:

```@example dofs
dhMixed = IGADofHandler(grid)
add!(dhMixed, :u, ip^2)
add!(dhMixed, :p, ip)
close!(dhMixed)

dhMixed.field_offsets, ndofs(dhMixed)
```

The fields are stored one after the other: all `:u` dofs first, then all `:p` dofs. Note that
every field uses the same spline basis here — a genuinely different basis per field (e.g. a
lower-order pressure space) requires a second geometry/basis object.

## Assembly

The assembly loop is a plain Ferrite loop. The only IGA-specific part is that `reinit!`
switches the active knot span of the interpolation:

```@example dofs
qr = QuadratureRule{RefQuadrilateral}(3)
cellvalues = CellValues(qr, ip^2, ip^2)

K = allocate_matrix(dh)
assembler = start_assemble(K)
ke = zeros(ndofs_per_cell(dh), ndofs_per_cell(dh))

for cell in CellIterator(dh)
    reinit!(cellvalues, cell)
    fill!(ke, 0.0)
    # ... element contributions, using shape_value/shape_gradient/getdetJdV as usual
    assemble!(assembler, celldofs(cell), ke)
end
nothing # hide
```

Note that the geometric interpolation passed to `CellValues` is the spline basis itself
(`ip^2` above, or `ip`), not a Lagrange interpolation: the geometry is the patch.

!!! note "IGAInterpolation is stateless"
    Older versions of this package stored the active knot span *in the interpolation
    object*, which made sharing one interpolation between tasks assembling different cells
    at the same time unsafe. `IGAInterpolation` no longer carries any such state — see
    [Parallel assembly](@ref) below.

## Parallel assembly

`IGAInterpolation` carries no per-cell state: `reinit!` remaps the quadrature points into
the active knot span before the (shared, read-only) spline basis ever sees them, instead of
storing "the current cell" anywhere. This makes it exactly as safe to share as a `Lagrange`
interpolation, so parallel assembly follows Ferrite's ordinary recipe — give every task its
own `CellCache`/`CellValues` for their private scratch buffers (e.g. with `TaskLocalValue`),
all built from the very same `ip`:

```julia
using Base.Threads: @threads
using TaskLocalValues: TaskLocalValue # as used in Ferrite's own threaded-assembly recipe

n = ndofs_per_cell(dh)
ccs = TaskLocalValue{CellCache}(() -> CellCache(dh))
cvs = TaskLocalValue{typeof(cellvalues)}(() -> copy(cellvalues))
kes = TaskLocalValue{Matrix{Float64}}(() -> zeros(n, n))

@threads for cellid in 1:getncells(dh)
    cc, cv, ke = ccs[], cvs[], kes[]
    reinit!(cc, cellid)
    reinit!(cv, cc)
    fill!(ke, 0.0)
    # ... element contributions, using shape_value/shape_gradient/getdetJdV as usual
    assemble!(assembler, celldofs(cc), ke) # a thread-safe assembler is required too
end
```

`copy(cellvalues)` here is only about giving each task its own `Nx`/`dNdx` scratch buffers —
the same requirement any Ferrite interpolation has for threaded assembly. There is nothing
IGA-specific to it any more: `ip` itself needs no copying at all, and reusing the very same
`ip` (or even the same `CellValues`, as long as its scratch buffers aren't written
concurrently) across tasks is fine.

## Evaluating a solution

Ferrite's `evaluate_at_grid_nodes` works and returns the values at the corners of the knot
spans (i.e. on `parameterSpaceGrid`). To evaluate anywhere else, use
[`interpolate`](@ref), which takes a *parametric* point:

```@example dofs
u = zeros(ndofs(dh))
interpolate(ip^2, u, [0.25, 0.5])
```

Passing the interpolation (not just the basis) is what makes it aware of the interleaved
layout; the scalar and the vector-valued case are handled by different methods. Use the
`offset` keyword to address a field that does not start at the beginning of the vector, e.g.
`interpolate(ip, u, x; offset = dhMixed.field_offsets[2])` for `:p` above.

## Limitations

- Only one `SubDofHandler` (one cell type, one patch) is supported.
- The handler relies on internal Ferrite functionality; see the
  [Developer documentation](../developer.md) for why, and for what a redesign would change.
