@testitem "evaluate_at_grid_nodes returns the field on the parameter-space grid" begin
    geometry = createBSplineSquare(1.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 2)
    grid = IGAGrid{2}(geometry)
    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))

    dh = IGADofHandler(grid)
    add!(dh, :u, ip^2)
    close!(dh)

    # An affine field in the control-point coordinates, cf. the interpolation tests
    u = zeros(ndofs(dh))
    for i in 1:getnnodes(grid)
        X = Ferrite.get_node_coordinate(grid, i)
        u[2 * (i - 1) + 1] = 3X[1] + 1
        u[2 * (i - 1) + 2] = -2X[2]
    end

    vals = evaluate_at_grid_nodes(dh, u, :u)
    paramGrid = parameterSpaceGrid(grid)

    # One value per knot-span corner — not one per control point — and none left as NaN
    @test length(vals) == getnnodes(paramGrid)
    @test length(vals) != getnnodes(grid)
    @test !any(v -> any(isnan, v), vals)

    for (nodeid, node) in enumerate(getnodes(paramGrid))
        ξ = node.x
        @test vals[nodeid][1] ≈ 3ξ[1] + 1 rtol = 1.0e-12
        @test vals[nodeid][2] ≈ -2ξ[2] rtol = 1.0e-12
    end
end

@testitem "L2 projection of quadrature point data" begin
    geometry = createBSplineSquare(1.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 2)
    grid = IGAGrid{2}(geometry)
    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    qr = QuadratureRule{RefQuadrilateral}(3)
    cellvalues = CellValues(qr, ip, ip)

    # The projector must be constructible without an explicit `qr_lhs`
    projector = L2Projector(ip, grid)
    @test projector isa Ferrite.AbstractProjector
    # A vector-valued interpolation projects onto its scalar base interpolation
    @test L2Projector(ip^2, grid) isa Ferrite.AbstractProjector

    # A linear function sampled at the quadrature points is reproduced exactly
    f(x) = 2x[1] - 3x[2] + 1
    qpdata = [zeros(getnquadpoints(qr)) for _ in 1:getncells(grid)]
    for cell in CellIterator(grid)
        reinit!(cellvalues, cell)
        for q_point in 1:getnquadpoints(cellvalues)
            x = spatial_coordinate(cellvalues, q_point, getcoordinates(cell))
            qpdata[cellid(cell)][q_point] = f(x)
        end
    end

    projected = project(projector, qpdata, qr)
    @test length(projected) == ndofs(projector.dh)

    vals = evaluate_at_grid_nodes(projector, projected)
    paramGrid = parameterSpaceGrid(grid)
    @test length(vals) == getnnodes(paramGrid)
    @test !any(isnan, vals)
    for (nodeid, node) in enumerate(getnodes(paramGrid))
        x = FerriteGismo.toPhysical(grid, Vector(node.x))
        @test vals[nodeid] ≈ f(x) rtol = 1.0e-10 atol = 1.0e-12   # f vanishes at some nodes
    end

    # Tensor-valued data works as well and can be written out
    σqp = [[SymmetricTensor{2, 2}((1.0, 0.0, 2.0)) for _ in 1:getnquadpoints(qr)] for _ in 1:getncells(grid)]
    σ = project(projector, σqp, qr)
    @test all(s -> s ≈ SymmetricTensor{2, 2}((1.0, 0.0, 2.0)), evaluate_at_grid_nodes(projector, σ))

    dh = IGADofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)
    mktempdir() do outdir
        VTKGridFile(joinpath(outdir, "projection"), dh) do vtk
            write_projection(vtk, projector, σ, "stress")
        end
        @test isfile(joinpath(outdir, "projection.vtu"))
    end
end

@testitem "The export mesh follows the knot spans" begin
    using LinearAlgebra

    geometry = createBSplineSquare(1.0)
    degreeElevate!(geometry, 1)
    # Deliberately non-uniform: the export mesh must follow the actual knot spans, not a
    # uniform subdivision of the parameter domain.
    insertKnot!(geometry, 0.1, 1)
    insertKnot!(geometry, 0.9, 1)
    insertKnot!(geometry, 0.5, 2)
    grid = IGAGrid{2}(geometry)

    @test breakpoints(grid, 1) ≈ [0.0, 0.1, 0.9, 1.0]
    @test breakpoints(grid, 2) ≈ [0.0, 0.5, 1.0]
    @test numElements(grid) == (3, 2)
    @test prod(numElements(grid)) == getncells(grid)

    eg = exportGrid(grid)
    @test eg isa Grid{2, Quadrilateral, Float64}
    @test getncells(eg) == getncells(grid)
    @test getnnodes(eg) == length(breakpoints(grid, 1)) * length(breakpoints(grid, 2))

    # The unit square is its own parameterization, so the node coordinates are the knots
    xs = sort(unique(round.([Ferrite.get_node_coordinate(eg, i)[1] for i in 1:getnnodes(eg)], digits = 12)))
    ys = sort(unique(round.([Ferrite.get_node_coordinate(eg, i)[2] for i in 1:getnnodes(eg)], digits = 12)))
    @test xs ≈ breakpoints(grid, 1)
    @test ys ≈ breakpoints(grid, 2)

    # The parametric mesh carries the same points, in the same order
    pg = parameterSpaceGrid(grid)
    @test getncells(pg) == getncells(eg)
    @test [n.x for n in getnodes(pg)] == exportPoints(grid)

    @test_throws ArgumentError breakpoints(grid, 3)
    @test_throws ArgumentError breakpoints(grid, 1; subdivision = 0)
end

@testitem "Subdivision resolves geometry and field" begin
    using LinearAlgebra

    r_in, r_out = 1.0, 2.0
    geometry = createNurbsQuarterAnnulus(r_in, r_out)
    uniformRefine!(geometry, 2)
    grid = IGAGrid{2}(geometry)
    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))

    dh = IGADofHandler(grid)
    add!(dh, :u, ip^2)
    close!(dh)

    # u ≡ the geometry map, so its values must equal the export node coordinates
    u = zeros(ndofs(dh))
    for i in 1:getnnodes(grid)
        X = Ferrite.get_node_coordinate(grid, i)
        u[2 * (i - 1) + 1], u[2 * (i - 1) + 2] = X[1], X[2]
    end

    for subdivision in (1, 3)
        eg = exportGrid(grid; subdivision)
        @test getncells(eg) == subdivision^2 * getncells(grid)

        vals = evaluateAtExportNodes(dh, u, :u; subdivision)
        @test length(vals) == getnnodes(eg)
        @test maximum(norm(vals[i] - Ferrite.get_node_coordinate(eg, i)) for i in 1:getnnodes(eg)) < 1.0e-12

        # Every export node lies on the exact annulus, also the subdivided ones
        radii = [norm(Ferrite.get_node_coordinate(eg, i)) for i in 1:getnnodes(eg)]
        @test all(r -> r_in - 1.0e-12 <= r <= r_out + 1.0e-12, radii)
    end

    # Subdividing samples the geometry between the knot-span corners
    eg3 = exportGrid(grid; subdivision = 3)
    @test any(r -> r_in + 0.05 < r < r_out - 0.05, [norm(Ferrite.get_node_coordinate(eg3, i)) for i in 1:getnnodes(eg3)])

    mktempdir() do outdir
        VTKGridFile(joinpath(outdir, "smooth"), dh; subdivision = 3) do vtk
            write_solution(vtk, dh, u; subdivision = 3)
        end
        @test isfile(joinpath(outdir, "smooth.vtu"))
    end
end

@testsnippet BezierEvaluation begin
    const FG = FerriteGismo

    # Rational Bézier representation of a cell evaluated at the local point `t` ∈ [0, 1]^rdim,
    # the way Paraview interpolates geometry and point data
    function evaluateBezier(degrees, coefs, weights, t)
        B = [FG._bernstein(degrees, idx, t) for idx in FG._bernsteinIndices(degrees)]
        return sum(weights .* B .* coefs) / sum(weights .* B)
    end

    toParameter(span, t) = [span.lower[d] + t[d] * (span.upper[d] - span.lower[d]) for d in eachindex(t)]

    # Largest deviation of the Bézier geometry and of the Bézier field of `u` from the spline
    function bezierErrors(grid, dh, u)
        rdim = Ferrite.get_reference_dimension(grid)
        basis = TinyGismo.basis(grid.geometry)
        weights = FG._rationalWeights(basis)
        fieldBasis = Ferrite.get_base_interpolation(FG._fieldInterpolation(dh, :u)).basis

        degrees = max.(FG._degrees(basis, Val(rdim)), FG._degrees(fieldBasis, Val(rdim)))
        ext = FG.BezierExtraction(degrees)
        geometry = [FG._bezierGeometry(ext, grid, basis, weights, span) for span in (c.span for c in grid.cells)]
        field = FG._bezierField(ext, grid, last.(geometry), fieldBasis, a -> u[a])

        nb = length(ext.points)
        errX = errU = 0.0
        for (e, span) in ((e, c.span) for (e, c) in pairs(grid.cells)), t in ([0.0, 0.0], [0.3, 0.8], [0.77, 0.21], [1.0, 1.0])
            t = t[1:rdim]
            xb, wb = geometry[e]
            ub = field[(nb * (e - 1) + 1):(nb * e)]
            ξ = toParameter(span, t)
            errX = max(errX, norm(evaluateBezier(degrees, xb, wb, t) - FG.toPhysical(grid, ξ)))
            errU = max(errU, abs(evaluateBezier(degrees, ub, wb, t) - interpolate(fieldBasis, u, ξ)))
        end
        return errX, errU
    end

    function scalarHandler(grid, ip)
        dh = IGADofHandler(grid)
        add!(dh, :u, ip)
        close!(dh)
        return dh
    end
end

@testitem "Bézier export reproduces geometry and fields exactly" setup = [BezierEvaluation] begin
    using LinearAlgebra

    square() = (g = createBSplineSquare(1.0); degreeElevate!(g); uniformRefine!(g, 2); g)

    @testset "B-spline patch" begin
        g = square()
        grid = IGAGrid{2}(g)
        dh = scalarHandler(grid, IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(g)))
        @test all(bezierErrors(grid, dh, sin.(1:ndofs(dh))) .< 1.0e-12)
    end

    @testset "NURBS patch with different degrees per direction" begin
        g = createNurbsQuarterAnnulus(1.0, 2.0)
        degreeElevate!(g)
        uniformRefine!(g, 2)
        grid = IGAGrid{2}(g)
        @test FG._degrees(TinyGismo.basis(g), Val(2)) == (2, 3)
        dh = scalarHandler(grid, IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(g)))
        @test all(bezierErrors(grid, dh, sin.(1:ndofs(dh))) .< 1.0e-12)
    end

    @testset "field of higher degree than the geometry" begin
        g = createBSplineSquare(1.0)
        uniformRefine!(g, 2)
        f = createBSplineSquare(1.0)
        degreeElevate!(f, 2)
        uniformRefine!(f, 2)
        grid = IGAGrid{2}(g)
        dh = scalarHandler(grid, IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(f)))
        @test all(bezierErrors(grid, dh, sin.(1:ndofs(dh))) .< 1.0e-12)
    end

    @testset "hierarchical patch" begin
        thb = THBSpline{2}(square())
        refineElements!(thb, RefinementBox(2, 1:2, 1:2))
        grid = IGAGrid{2}(thb)
        dh = IGADofHandler(grid)
        for (cells, ip) in hierarchicalSubdomains(grid)
            add!(SubDofHandler(dh, cells), :u, ip)
        end
        close!(dh)
        @test all(bezierErrors(grid, dh, sin.(1:ndofs(dh))) .< 1.0e-12)
    end

    @testset "NURBS curve" begin
        c = createNurbsCircle(1.0)
        grid = IGAGrid{2, 1}(c)
        dh = scalarHandler(grid, IGAInterpolation{RefLine}(TinyGismo.basis(c)))
        @test all(bezierErrors(grid, dh, sin.(1:ndofs(dh))) .< 1.0e-12)
    end
end

@testitem "Bézier node ordering follows VTK" begin
    using FerriteGismo: _vtkOrdering
    # vertices, then the interior nodes, of a curve
    @test _vtkOrdering((3,)) == [1, 4, 2, 3]
    # vertices counterclockwise, edges in the direction of increasing parameter, then the interior
    @test _vtkOrdering((2, 2)) == [1, 3, 9, 7, 2, 6, 8, 4, 5]
    @test _vtkOrdering((3, 1)) == [1, 4, 8, 5, 2, 3, 6, 7]
end

@testitem "VTKBezierFile" setup = [BezierEvaluation] begin
    g = createNurbsQuarterAnnulus(1.0, 2.0)
    degreeElevate!(g)
    uniformRefine!(g, 2)
    grid = IGAGrid{2}(g)
    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(g))
    dh = IGADofHandler(grid)
    add!(dh, :u, ip^2)
    close!(dh)
    u = rand(ndofs(dh))

    proj = L2Projector(ip, grid)
    qr = QuadratureRule{RefQuadrilateral}(3)
    projected = project(proj, [ones(getnquadpoints(qr)) for _ in 1:getncells(grid)], qr)

    mktempdir() do outdir
        pvd = FG.WriteVTK.paraview_collection(joinpath(outdir, "series"))
        VTKBezierFile(joinpath(outdir, "annulus"), dh) do vtk
            write_solution(vtk, dh, u)
            write_projection(vtk, proj, projected, "projected")
            write_cell_data(vtk, collect(1.0:getncells(grid)), "cellid")
            pvd[0.0] = vtk
        end
        FG.WriteVTK.vtk_save(pvd)

        vtu = read(joinpath(outdir, "annulus.vtu"), String)
        @test isfile(joinpath(outdir, "series.pvd"))
        @test occursin("RationalWeights=\"RationalWeights\"", vtu)
        @test occursin("HigherOrderDegrees=\"HigherOrderDegrees\"", vtu)
        @test all(name -> occursin("Name=\"$name\"", vtu), ("u", "projected", "cellid"))

        # A B-spline field on a NURBS geometry has another weight function than the cells
        b = createBSplineSquare(1.0)
        degreeElevate!(b)
        uniformRefine!(b, 2)
        VTKBezierFile(joinpath(outdir, "invalid"), dh) do vtk
            @test_throws ArgumentError FG._bezierField(vtk, TinyGismo.basis(b), a -> 1.0)
        end
    end
end
