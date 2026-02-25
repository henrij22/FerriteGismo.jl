@testitem "CellValues Test" begin
    geometry = createBSplineSquare(2.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 1)
    grid = IGAGrid{2}(geometry)

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    qr = QuadratureRule{RefQuadrilateral}(1)

    cv = CellValues(qr, ip, ip; update_hessians = true)

    dh = IGADofHandler(grid)
    add!(dh, :w, ip)
    close!(dh)

    cc = CellCache(dh)
    reinit!(cc, 1)
    reinit!(cv, cc)

    @test getnbasefunctions(cv) == 9
    @test Ferrite.function_difforder(cv) == 2
    @test shape_value(cv, 1, 1) == 0.0625
    @test shape_gradient(cv, 1, 1) == [-0.25, -0.25]
    @test shape_hessian(cv, 1, 1) == [0.5 1.0; 1.0 0.5]

    @test typeof(shape_value(cv, 1, 1)) == Ferrite.shape_value_type(cv)
    @test typeof(shape_gradient(cv, 1, 1)) == Ferrite.shape_gradient_type(cv)
    @test typeof(shape_hessian(cv, 1, 1)) == Ferrite.shape_hessian_type(cv)

    @test typeof(Ferrite.function_interpolation(cv)) == typeof(ip)
    @test getnquadpoints(cv) == 1
    @test getdetJdV(cv, 1) ≈ 1.0
end
