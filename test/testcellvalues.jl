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

@testitem "Hessians on embedded elements (curve in 2D)" begin
    # u = x₁ on the unit circle: ∇ₛ∇ₛu = -(t₁ n + n₁ t) ⊗ t
    geometry = createNurbsCircle(1.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 1)
    grid = IGAGrid{2, 1}(geometry)
    ip = IGAInterpolation{RefLine}(TinyGismo.basis(geometry))
    cv = CellValues(QuadratureRule{RefLine}(4), ip, ip^2; update_hessians = true)

    dh = IGADofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)

    u = [node.x[1] for node in grid.nodes]
    cc = CellCache(dh)
    for cell in 1:getncells(grid)
        reinit!(cc, cell)
        reinit!(cv, cc)
        u_cell = u[celldofs(cc)]
        for q in 1:getnquadpoints(cv)
            n = spatial_coordinate(cv, q, getcoordinates(cc))
            @test norm(n) ≈ 1.0
            t = Vec{2}((-n[2], n[1]))
            expected = -(t[1] * n + n[1] * t) ⊗ t
            @test Ferrite.function_hessian(cv, q, u_cell) ≈ expected atol = 1.0e-10
        end
    end
end
