@testitem "prescribeEdge! on a vector field" begin
    using Ferrite

    geometry = createBSplineSquare(1.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 2)
    grid = IGAGrid{2}(geometry)

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    dh = IGADofHandler(grid)
    add!(dh, :u, ip^2)
    close!(dh)

    ch = ConstraintHandler(dh)
    prescribeEdge!(dh, ch, :right, :u, (x, t) -> 0.1t; components = 1)
    fixEdge!(dh, ch, :left, :u)
    close!(ch)

    left = edgeControlPoints(dh, :left, :u)
    right = edgeControlPoints(dh, :right, :u)
    @test length(right) == length(left) > 1
    # One component prescribed on the right edge, both on the left one
    @test length(ch.prescribed_dofs) == length(right) + 2 * length(left)

    u = zeros(ndofs(dh))
    Ferrite.update!(ch, 2.0)
    apply!(u, ch)
    @test all(u[2 * (p - 1) + 1] ≈ 0.2 for p in right)
    @test all(iszero(u[2 * (p - 1) + c]) for p in left, c in 1:2)
    # The unconstrained component of the right edge is untouched
    @test all(iszero(u[2 * (p - 1) + 2]) for p in right)

    # A number instead of a function, and a condition without `components`
    ch2 = ConstraintHandler(dh)
    prescribeEdge!(dh, ch2, :top, :u, 0.5)
    close!(ch2)
    u2 = zeros(ndofs(dh))
    apply!(u2, ch2)
    top = edgeControlPoints(dh, :top, :u)
    @test length(ch2.prescribed_dofs) == 2 * length(top)
    @test all(u2[2 * (p - 1) + c] ≈ 0.5 for p in top, c in 1:2)
end

@testitem "Prescribed displacement of a whole edge is met exactly" begin
    using Ferrite, SparseArrays

    # Uniaxial stretch of the unit square: prescribing the same value on every control
    # point of an edge reproduces a constant edge displacement exactly (partition of unity).
    geometry = createBSplineSquare(1.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 2)
    grid = IGAGrid{2}(geometry)

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    qr = QuadratureRule{RefQuadrilateral}(3)
    cellvalues = CellValues(qr, ip^2, ip^2)

    dh = IGADofHandler(grid)
    add!(dh, :u, ip^2)
    close!(dh)

    ch = ConstraintHandler(dh)
    fixEdge!(dh, ch, :left, :u; components = 1)
    fixEdge!(dh, ch, :bottom, :u; components = 2)
    prescribeEdge!(dh, ch, :right, :u, (x, t) -> 0.1t; components = 1)
    close!(ch)
    Ferrite.update!(ch, 1.0)

    Emod, ν = 200.0e3, 0.3
    Gmod = Emod / (2 * (1 + ν))
    Kmod = Emod / (3 * (1 - 2ν))
    C = gradient(ϵ -> 2 * Gmod * dev(ϵ) + 3 * Kmod * vol(ϵ), zero(SymmetricTensor{2, 2}))

    K = allocate_matrix(dh)
    assembler = start_assemble(K)
    ke = zeros(getnbasefunctions(cellvalues), getnbasefunctions(cellvalues))
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        fill!(ke, 0.0)
        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)
            for i in 1:getnbasefunctions(cellvalues), j in 1:getnbasefunctions(cellvalues)
                ke[i, j] += (shape_gradient(cellvalues, q_point, i) ⊡ C ⊡ shape_symmetric_gradient(cellvalues, q_point, j)) * dΩ
            end
        end
        assemble!(assembler, celldofs(cell), ke)
    end

    f = zeros(ndofs(dh))
    apply!(K, f, ch)
    u = K \ f
    apply!(u, ch)

    # u_x is linear in x and reaches the prescribed value on the whole right edge
    for ξ in (0.0, 0.25, 0.5, 0.75, 1.0)
        @test interpolate(ip^2, u, [1.0, ξ])[1] ≈ 0.1 rtol = 1.0e-10
        @test interpolate(ip^2, u, [0.0, ξ])[1] ≈ 0.0 atol = 1.0e-12
        @test interpolate(ip^2, u, [0.5, ξ])[1] ≈ 0.05 rtol = 1.0e-8
    end
end

@testitem "edgeControlPoints picks up the control points of each side" begin
    using Ferrite

    geometry = createBSplineSquare(1.0)
    degreeElevate!(geometry, 1)      # quadratic
    uniformRefine!(geometry, 2)      # 3 x 3 elements -> 5 x 5 control points
    grid = IGAGrid{2}(geometry)

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    dh = IGADofHandler(grid)
    add!(dh, :u, ip^2)
    close!(dh)

    sides = (:left, :right, :bottom, :top)
    cps = Dict(side => edgeControlPoints(dh, side, :u) for side in sides)

    n = isqrt(getnnodes(grid))
    @test getnnodes(grid) == n^2
    for side in sides
        @test length(cps[side]) == n
        @test all(1 .<= collect(cps[side]) .<= getnnodes(grid))
    end
    # Opposite sides are disjoint, adjacent ones share exactly the corner control point
    @test isempty(intersect(cps[:left], cps[:right]))
    @test isempty(intersect(cps[:bottom], cps[:top]))
    for a in (:left, :right), b in (:bottom, :top)
        @test length(intersect(cps[a], cps[b])) == 1
    end
    # The union of all sides is the boundary layer of the control net
    @test length(union(values(cps)...)) == 4n - 4

    # The control points really lie on the corresponding side of the patch
    for (side, coord, value) in ((:left, 1, 0.0), (:right, 1, 1.0), (:bottom, 2, 0.0), (:top, 2, 1.0))
        @test all(Ferrite.get_node_coordinate(grid, p)[coord] ≈ value for p in cps[side])
    end

    # The G+Smo boundary index may be given instead of the symbol
    @test edgeControlPoints(dh, 1, :u) == cps[:left]
end

@testitem "prescribeEdge! follows the time / load factor" begin
    using Ferrite

    geometry = createBSplineSquare(1.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 2)
    grid = IGAGrid{2}(geometry)

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    dh = IGADofHandler(grid)
    add!(dh, :u, ip^2)
    close!(dh)

    ch = ConstraintHandler(dh)
    prescribeEdge!(dh, ch, :left, :u, (x, t) -> -0.5t; components = 1)
    close!(ch)

    right = edgeControlPoints(dh, :left, :u)
    dofs = [2 * (p - 1) + 1 for p in right]
    for t in (0.0, 1.0, 3.5)
        Ferrite.update!(ch, t)
        u = zeros(ndofs(dh))
        apply!(u, ch)
        @test all(u[d] ≈ -0.5t for d in dofs)

        # apply_zero! must not write the inhomogeneity
        u0 = ones(ndofs(dh))
        apply_zero!(u0, ch)
        @test all(iszero(u0[d]) for d in dofs)
    end

    # A one-argument f(x) is constant in time, as for a plain Ferrite Dirichlet
    ch2 = ConstraintHandler(dh)
    prescribeEdge!(dh, ch2, :left, :u, x -> 0.25; components = 2)
    close!(ch2)
    for t in (0.0, 2.0)
        Ferrite.update!(ch2, t)
        u = zeros(ndofs(dh))
        apply!(u, ch2)
        @test all(u[2 * (p - 1) + 2] ≈ 0.25 for p in right)
    end
end

@testitem "prescribeEdge! evaluates f at the control-point coordinates" begin
    using Ferrite

    geometry = createBSplineSquare(1.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 2)
    grid = IGAGrid{2}(geometry)

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    dh = IGADofHandler(grid)
    add!(dh, :u, ip^2)
    close!(dh)

    f(x, t) = t * x[2]^2
    ch = ConstraintHandler(dh)
    prescribeEdge!(dh, ch, :right, :u, f; components = 1)
    close!(ch)
    Ferrite.update!(ch, 2.0)

    u = zeros(ndofs(dh))
    apply!(u, ch)
    for p in edgeControlPoints(dh, :right, :u)
        @test u[2 * (p - 1) + 1] ≈ f(Ferrite.get_node_coordinate(grid, p), 2.0)
    end
end

@testitem "A zero prescribeEdge! constrains the same dofs as fixEdge!" begin
    using Ferrite

    geometry = createBSplineSquare(1.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 1)
    grid = IGAGrid{2}(geometry)

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    dh = IGADofHandler(grid)
    add!(dh, :u, ip^2)
    close!(dh)

    for components in (nothing, 1, [1, 2])
        chFix = ConstraintHandler(dh)
        fixEdge!(dh, chFix, :top, :u; components)
        close!(chFix)

        chPrescribe = ConstraintHandler(dh)
        prescribeEdge!(dh, chPrescribe, :top, :u, 0.0; components)
        close!(chPrescribe)

        @test chPrescribe.prescribed_dofs == chFix.prescribed_dofs
        @test all(iszero, chPrescribe.inhomogeneities)
    end
end

@testitem "prescribeEdge! on a scalar field" begin
    using Ferrite

    geometry = createBSplineSquare(1.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 1)
    grid = IGAGrid{2}(geometry)

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    dh = IGADofHandler(grid)
    add!(dh, :T, ip)
    close!(dh)

    ch = ConstraintHandler(dh)
    prescribeEdge!(dh, ch, :bottom, :T, 100.0)
    close!(ch)

    bottom = edgeControlPoints(dh, :bottom, :T)
    @test ch.prescribed_dofs == sort(collect(bottom))   # one dof per control point
    u = zeros(ndofs(dh))
    apply!(u, ch)
    @test all(u[p] ≈ 100.0 for p in bottom)
    # A constant prescribed on a whole edge is reproduced exactly in between
    @test interpolate(ip, u, [0.5, 0.0]) ≈ 100.0
end

@testitem "Invalid IGA Dirichlet conditions are rejected" begin
    using Ferrite

    geometry = createBSplineSquare(1.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 1)
    grid = IGAGrid{2}(geometry)

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    dh = IGADofHandler(grid)
    add!(dh, :u, ip^2)
    close!(dh)

    # Facet sets carry no dof information on an IGA grid, only control points do
    ch = ConstraintHandler(dh)
    @test_throws ArgumentError add!(ch, Dirichlet(:u, getfacetset(grid, "left"), (x, t) -> 0.0, 1))

    # Components outside the field
    @test_throws ArgumentError prescribeEdge!(dh, ch, :left, :u, 0.0; components = 3)

    # Unknown side and unknown field
    @test_throws ArgumentError prescribeEdge!(dh, ch, :front, :u, 0.0)
    @test_throws ErrorException prescribeEdge!(dh, ch, :left, :T, 0.0)

    # Nothing may be added after closing
    close!(ch)
    @test_throws ArgumentError prescribeEdge!(dh, ch, :right, :u, 0.0)
end
