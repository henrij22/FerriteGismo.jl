@testitem "Grid Test" begin
    geometry = createBSplineSquare(2.0)
    grid = IGAGrid{2}(geometry)

    @test length(getcells(grid)) == 1
    @test length(first(getcells(grid)).nodes) == 4
    @test first(getcells(grid)).nodes == collect(1:4)

    degreeElevate!(geometry, 1)
    grid = IGAGrid{2}(geometry)

    @test length(getcells(grid)) == 1
    @test length(first(getcells(grid)).nodes) == 9
    @test first(getcells(grid)).nodes == collect(1:9)
    @test length(getnodes(grid)) == 9

    uniformRefine!(geometry, 1)
    grid = IGAGrid{2}(geometry)

    @test length(getcells(grid)) == 4
    @test length(first(getcells(grid)).nodes) == 9
    @test length(getnodes(grid)) == 16

    # Transformations
    @test FerriteGismo.toPhysical(grid, [1.0, 1.0]) ≈ Vec{2}((2.0, 2.0))
    @test FerriteGismo.ref_to_param([0.0, 0.0], first(grid.knotSpans)) ≈ [0.25, 0.25]

    # ParameterSpaceGrid
    uniformgrid = parameterSpaceGrid(grid)
    @test length(getcells(grid)) == 4
end
