# Applying forces

Neumann boundary conditions and body loads need no IGA-specific API: they are integrals over
facets and cells, and both `FacetValues` and `CellValues` support spline bases. The only IGA
detail is that the facet sets of a patch come from its parametric sides.

## Traction on a boundary

```@example forces
using FerriteGismo

geometry = createBSplineSquare(1.0)
degreeElevate!(geometry, 1)
uniformRefine!(geometry, 2)
grid = IGAGrid{2}(geometry)

ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
fqr = FacetQuadratureRule{RefQuadrilateral}(3)
facetvalues = FacetValues(fqr, ip^2, ip^2)

dh = IGADofHandler(grid)
add!(dh, :u, ip^2)
close!(dh)

traction = Vec{2}((0.0, -1.0e3))

f = zeros(ndofs(dh))
fe = zeros(getnbasefunctions(facetvalues))
for facet in FacetIterator(dh, getfacetset(grid, "right"))
    reinit!(facetvalues, facet)
    fill!(fe, 0.0)
    for q_point in 1:getnquadpoints(facetvalues)
        dΓ = getdetJdV(facetvalues, q_point)
        for i in 1:getnbasefunctions(facetvalues)
            fe[i] += (traction ⋅ shape_value(facetvalues, q_point, i)) * dΓ
        end
    end
    assemble!(f, celldofs(facet), fe)
end

sum(f[2:2:end])   # total applied force = traction * edge length
```

That is the standard Ferrite loop, unchanged. `getdetJdV` already contains the mapping from
the reference facet through the parametric facet to the physical boundary, and
`getnormal(facetvalues, q_point)` gives the outward normal — on a curved patch (a NURBS
annulus, say) it follows the exact geometry, not a polygonal approximation.

For a load that varies over the boundary, evaluate the physical coordinate at the quadrature
point:

```julia
x = spatial_coordinate(facetvalues, q_point, getcoordinates(facet))
tₚ = traction(x)
```

Note the asymmetry with Dirichlet conditions: a traction is integrated at quadrature points,
so a spatially varying Neumann condition is captured consistently, whereas a varying
Dirichlet condition is only interpolated at the control points (see
[Boundary conditions](boundary_conditions.md)).

## Body loads

The same holds for volumetric loads, integrated with `CellValues` over the cells:

```@example forces
qr = QuadratureRule{RefQuadrilateral}(3)
cellvalues = CellValues(qr, ip^2, ip^2)
b = Vec{2}((0.0, -9.81))

fb = zeros(ndofs(dh))
fe = zeros(ndofs_per_cell(dh))
for cell in CellIterator(dh)
    reinit!(cellvalues, cell)
    fill!(fe, 0.0)
    for q_point in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:getnbasefunctions(cellvalues)
            fe[i] += (b ⋅ shape_value(cellvalues, q_point, i)) * dΩ
        end
    end
    assemble!(fb, celldofs(cell), fe)
end
nothing # hide
```

## Point forces

There is no equivalent of a nodal force: control points are not points of the geometry, so
writing a value into a control-point dof does **not** apply that force at that location. A
point load has to be applied as the corresponding load functional, i.e. by evaluating the
basis at the parametric point ``\xi`` of the load and distributing it over the active
functions:

```@example forces
using FerriteGismo.TinyGismo: active!, eval!

basis = TinyGismo.basis(geometry)
ξ = [0.5, 1.0]                   # parametric position of the load
P = Vec{2}((0.0, -5.0))

actives = gsMatrix{Int32}()
values = gsMatrix()
active!(basis, ξ, actives)
eval!(basis, ξ, values)

fp = zeros(ndofs(dh))
for (k, p) in enumerate(toVector(actives))
    N = toVector(values)[k]
    for c in 1:2
        fp[2 * (p - 1) + c] += N * P[c]
    end
end
sum(fp[2:2:end])                 # equals P[2]
```

The load ends up spread over the ``(p+1)^2`` control points whose functions are non-zero at
``\xi``, which is the correct discrete representation of a point load in the spline space.

## Quadrature order

Spline elements are not Lagrange elements: a knot span carries basis functions of degree `p`,
so a `QuadratureRule` of order `p + 1` is the usual choice, and the cost per element is
higher than for a Lagrange element of the same degree. Reduced or patch-wise quadrature
rules are not provided.
