@testitem "Interpolation Test" begin
    geometry = createBSplineSquare(2.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 1)
    grid = IGAGrid{2}(geometry)

    nip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    @test getnbasefunctions(nip) == 9 # Per cell

    ipVec = nip^1
    @test Ferrite.n_components(ipVec) == 1
    @test Ferrite.get_base_interpolation(ipVec) == nip
    @test getnbasefunctions(ipVec) == 9

    ipVec = nip^3
    @test Ferrite.n_components(ipVec) == 3
    @test Ferrite.get_base_interpolation(ipVec) == nip
    @test getnbasefunctions(ipVec) == 27
end

@testitem "interpolate of a scalar IGA field" begin
    using Ferrite

    geometry = createBSplineSquare(1.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 2)
    grid = IGAGrid{2}(geometry)
    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    basis = TinyGismo.basis(geometry)

    # The interpolation method must agree with the `gsBasis` one it wraps
    u = collect(1.0:TinyGismo.size(basis))
    for x in ([0.0, 0.0], [0.3, 0.7], [1.0, 1.0])
        @test interpolate(ip, u, x) == interpolate(basis, u, x)
    end

    # Partition of unity: a constant coefficient vector is reproduced everywhere
    uc = fill(2.5, Int(TinyGismo.size(basis)))
    for x in ([0.0, 0.5], [0.25, 0.25], [1.0, 0.75])
        @test interpolate(ip, uc, x) ≈ 2.5
    end

    # A field given by the control-point coordinates is the (linear) geometry itself
    ux = [Ferrite.get_node_coordinate(grid, i)[1] for i in 1:getnnodes(grid)]
    for ξ in ([0.0, 0.0], [0.3, 0.7], [0.6, 0.1], [1.0, 1.0])
        @test interpolate(ip, ux, ξ) ≈ ξ[1] rtol = 1.0e-12
    end
end

@testitem "interpolate of a vector-valued IGA field" begin
    using Ferrite

    geometry = createBSplineSquare(1.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 2)
    grid = IGAGrid{2}(geometry)
    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))

    dh = IGADofHandler(grid)
    add!(dh, :u, ip^2)
    close!(dh)

    # Interleaved dofs u1x, u1y, u2x, u2y, ...: an affine field in the control-point
    # coordinates is reproduced exactly (the control points are the Greville abscissae
    # of this uniform patch).
    u = zeros(ndofs(dh))
    for i in 1:getnnodes(grid)
        X = Ferrite.get_node_coordinate(grid, i)
        u[2 * (i - 1) + 1] = 3X[1] + 1
        u[2 * (i - 1) + 2] = -2X[2]
    end

    for ξ in ([0.0, 0.0], [0.25, 0.75], [0.5, 0.5], [1.0, 1.0])
        val = interpolate(ip^2, u, ξ)
        @test val isa Vec{2, Float64}
        @test val[1] ≈ 3ξ[1] + 1 rtol = 1.0e-12
        @test val[2] ≈ -2ξ[2] rtol = 1.0e-12
        # A `Vec` may be passed instead of a plain vector
        @test interpolate(ip^2, u, Vec{2}((ξ[1], ξ[2]))) ≈ val
    end

    # The components are picked apart correctly, i.e. it is not the scalar method in disguise
    @test interpolate(ip^2, u, [0.5, 0.5])[1] != interpolate(ip^2, u, [0.5, 0.5])[2]

    # `offset` addresses a field further back in the dof vector
    ushifted = vcat(zeros(7), u)
    @test interpolate(ip^2, ushifted, [0.25, 0.75]; offset = 7) ≈ interpolate(ip^2, u, [0.25, 0.75])

    # Three components work as well
    dh3 = IGADofHandler(grid)
    add!(dh3, :u, ip^3)
    close!(dh3)
    u3 = zeros(ndofs(dh3))
    for i in 1:getnnodes(grid)
        X = Ferrite.get_node_coordinate(grid, i)
        u3[3 * (i - 1) + 1] = X[1]
        u3[3 * (i - 1) + 2] = X[2]
        u3[3 * (i - 1) + 3] = 1.0
    end
    val3 = interpolate(ip^3, u3, [0.2, 0.4])
    @test val3 isa Vec{3, Float64}
    @test val3 ≈ Vec{3}((0.2, 0.4, 1.0)) rtol = 1.0e-12
end
