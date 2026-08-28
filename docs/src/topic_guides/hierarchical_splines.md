# Hierarchical splines

A tensor-product spline space cannot be refined locally. Inserting a knot inserts it across
the whole patch, so resolving a corner singularity costs elements everywhere along two
lines. Hierarchical splines break that coupling: a coarse tensor basis is the first *level*,
and finer levels are introduced only over the regions that need them.

FerriteGismo supports both flavours that G+Smo provides, through
[TinyGismo.jl](https://github.com/henrij22/TinyGismo.jl):

| | Partition of unity after refinement | |
|:--|:--|:--|
| `THBSplineBasis` / `THBSpline` | yes | truncated, the usual choice |
| `HBSplineBasis` / `HBSpline` | no | untruncated |

Everything on this page applies to both. The one place they differ is discussed under
[Truncation](@ref).

## Building a locally refined patch

Start from an ordinary tensor patch and lift it. The lift changes nothing geometrically — it
only makes the patch refinable element by element:

```@example hier
using FerriteGismo

geometry = createBSplineSquare(1.0)
degreeElevate!(geometry)             # quadratic
uniformRefine!(geometry, 2)          # 3x3 coarse elements

patch = THBSpline{2}(geometry)
numLevels(TinyGismo.basis(patch))
```

Refinement is addressed with a TinyGismo's `RefinementBox`, which
names the level a region is refined **to** and indexes its corners on *that level's* grid.
Each level halves the cells of the one before it, so the first coarse cell is cells `1:2` on
level 2:

```@example hier
refineElements!(patch, RefinementBox(2, 1:2, 1:2))   # coarse cell (1,1) -> level 2
refineElements!(patch, RefinementBox(3, 1:2, 1:2))   # its own first quarter -> level 3
numLevels(TinyGismo.basis(patch))
```

The grid is built exactly as for a tensor patch:

```@example hier
grid = IGAGrid{2}(patch)
getncells(grid), isHierarchical(grid)
```

`getLevelAtPoint` confirms where the refinement landed:

```@example hier
[getLevelAtPoint(TinyGismo.basis(patch), p) for p in ([0.05, 0.05], [0.2, 0.2], [0.8, 0.8])]
```

Drawing [`parameterSpaceGrid`](@ref) shows the graded mesh. On a hierarchical grid this mesh
is built element by element, so it is the real element layout and not a lattice over it:

```@example hier
using FerriteTikz
TikzPictures.tikzUseTectonic(true) # hide

tikzgrid(parameterSpaceGrid(grid); cellcolor = "blue!12", picturescale = 6.0)
```

Three levels are visible in the lower-left corner, and the rest of the patch is untouched.

## Why one interpolation is not enough

On a tensor patch every element sees the same number of basis functions, ``(p+1)^{d}``. That
is what lets a single [`IGAInterpolation`](@ref) describe the whole patch: Ferrite's
`getnbasefunctions` is one number per interpolation, and it sizes every `CellValues` buffer.

A hierarchical patch breaks that. An element next to a level transition sees the fine
functions on it *and* the coarse functions that overlap it:

```@example hier
sort(unique(length(cell.nodes) for cell in grid.cells))
```

So `IGAInterpolation` refuses a hierarchical basis rather than silently picking one of these
counts:

```@example hier
try
    IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(patch))
catch e
    println(e.msg)
end
```

## Grouping elements into subdomains

The fix is the mechanism Ferrite already has for "different interpolations on different
cells": a `SubDofHandler` per group. [`hierarchicalSubdomains`](@ref) partitions the cells by
their active count and hands back one interpolation per group:

```@example hier
groups = hierarchicalSubdomains(grid)
[(getnbasefunctions(ip), length(cells)) for (cells, ip) in groups]
```

Colouring the mesh by group shows what the partition means geometrically — the interior of
each level is uniform, and the extra functions live in a band along the level transitions:

```@example hier
mesh = parameterSpaceGrid(grid)
for (i, (cells, _)) in enumerate(groups)
    addcellset!(mesh, "group$i", Set(cells))
end

palette = ["blue!12", "orange!25", "teal!25", "purple!20", "red!20", "green!20"]
tikzgrid(
    mesh;
    cellsetcolors = ["group$i" => palette[mod1(i, length(palette))] for i in eachindex(groups)],
    picturescale = 6.0,
)
```

Each group gets its own `SubDofHandler`:

```@example hier
dh = IGADofHandler(grid)
for (cells, ip) in groups
    sdh = SubDofHandler(dh, cells)
    add!(sdh, :u, ip)
end
close!(dh)

ndofs(dh), Int(size(TinyGismo.basis(patch)))
```

The dof count equals the number of control points of the patch. That is the important
property: the subdomains partition the *elements*, not the basis, so dofs are numbered
globally over the control points and a control point active in two groups is the same dof in
both. Nothing about the subdivision leaks into the linear system.

## Assembling

The assembly loop gains exactly one level of nesting: iterate `dh.subdofhandlers`, and build
a `CellValues` per subdomain, since each has its own number of shape functions.

```julia
for sdh in dh.subdofhandlers
    ip = only(sdh.field_interpolations)
    cellvalues = CellValues(qr, ip, ip^2)
    n = getnbasefunctions(ip)
    # ... element matrices of size n, assembled with celldofs(cell) as usual
    for cell in CellIterator(sdh)
        reinit!(cellvalues, cell)
        # ...
    end
end
```

On a tensor patch there is a single subdomain and this reduces to the usual loop, so the
same code covers both. The
[Locally refined heat equation](../tutorials/hierarchical_refinement.md) tutorial works this
through end to end, including boundary conditions and the solve.

Boundary conditions need no special handling at all: constraints are expressed on control
points, which belong to the patch rather than to any subdomain, so [`fixEdge!`](@ref) and
[`prescribeEdge!`](@ref) behave exactly as on a tensor patch.

## Visualization

[`exportGrid`](@ref) and [`evaluateAtExportNodes`](@ref) work unchanged, and on a
hierarchical grid they too are built element by element:

```@example hier
getncells(exportGrid(grid)), getncells(grid)
```

Nodes are duplicated along element boundaries as a result — each element carries its own
corners. That is harmless for drawing and for pointwise field values, and it is what lets
the graded mesh be drawn exactly. `subdivision` samples each element several times per
direction, which is what you want for a curved geometry or a high-order field:

```@example hier
getncells(exportGrid(grid; subdivision = 3))
```

Two functions are *not* available, because they presuppose a tensor-product element lattice
that a hierarchical patch does not have. Both raise rather than answer wrongly:

- `breakpoints(grid, dir)` — there is no single set of breakpoints per direction.
- `numElementsPerDirection(grid, dir)` — the elements do not form a grid. Use
  `getncells(grid)` for the total.

## Truncation

`THBSplineBasis` truncates the coarse functions against the finer levels; `HBSplineBasis`
does not. Both span the same space and give the same number of dofs, but only the truncated
one remains a partition of unity over refined regions:

```@example hier
pou(patch, ξ) = sum(toMatrix(_eval(TinyGismo.basis(patch), ξ)))

truncated = THBSpline{2}(geometry)
untruncated = HBSpline{2}(geometry)
refineElements!(truncated, RefinementBox(2, 1:2, 1:2))
refineElements!(untruncated, RefinementBox(2, 1:2, 1:2))

inside = [0.05, 0.05]     # inside the refined region
pou(truncated, inside), pou(untruncated, inside)
```

Partition of unity is what makes the coefficients behave like control points and keeps the
system well-conditioned, so `THBSpline` is the default choice. `HBSpline` goes through the
identical FerriteGismo pipeline if you want the untruncated functions.

!!! note "Admissible refinement is not enforced"
    `refineElements!` refines exactly the elements it is given. It does not enforce a bound
    on how many levels may meet at an element, which analysis-suitability results rely on.
    G+Smo's `gsHBox` machinery for admissible refinement is not wrapped by TinyGismo yet, so
    for now grade your refinement yourself — refine each level over a region contained in the
    previous one, as the example on this page does.
