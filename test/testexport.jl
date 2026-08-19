@testitem "evaluate_at_grid_nodes returns the field on the parameter-space grid" begin
    using Ferrite

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
    using Ferrite

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
    using Ferrite, LinearAlgebra

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
    using Ferrite, LinearAlgebra

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
