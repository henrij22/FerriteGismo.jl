@testitem "FacetValues on the unit square" begin
    using Ferrite

    geometry = createBSplineSquare(1.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 2)
    grid = IGAGrid{2}(geometry)

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    fqr = FacetQuadratureRule{RefQuadrilateral}(3)
    fv = FacetValues(fqr, ip, ip)

    dh = IGADofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)

    expected_normals = Dict(
        "left" => Vec((-1.0, 0.0)), "right" => Vec((1.0, 0.0)),
        "bottom" => Vec((0.0, -1.0)), "top" => Vec((0.0, 1.0)),
    )

    for side in ("left", "right", "bottom", "top")
        facetset = getfacetset(grid, side)
        @test length(facetset) > 0
        len = 0.0
        for facet in FacetIterator(dh, facetset)
            reinit!(fv, facet)
            for q_point in 1:getnquadpoints(fv)
                len += getdetJdV(fv, q_point)
                @test getnormal(fv, q_point) ≈ expected_normals[side]
                # Partition of unity on the facet
                @test sum(shape_value(fv, q_point, i) for i in 1:getnbasefunctions(fv)) ≈ 1.0
            end
        end
        # Each edge of the unit square has length 1
        @test len ≈ 1.0
    end

    @test_throws ArgumentError getfacetset(grid, "front")
end

@testitem "FacetValues on a curved boundary (quarter annulus)" begin
    using Ferrite, LinearAlgebra

    r_in, r_out = 1.0, 2.0
    geometry = createNurbsQuarterAnnulus(r_in, r_out)
    uniformRefine!(geometry, 3)
    grid = IGAGrid{2}(geometry)

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    fqr = FacetQuadratureRule{RefQuadrilateral}(4)
    fv = FacetValues(fqr, ip, ip)

    dh = IGADofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)

    # The parametric "left"/"right" sides map to the inner/outer circular arc.
    expected_lengths = Dict(
        "left" => r_in * π / 2, "right" => r_out * π / 2,
        "bottom" => r_out - r_in, "top" => r_out - r_in,
    )

    for side in ("left", "right", "bottom", "top")
        len = 0.0
        for facet in FacetIterator(dh, getfacetset(grid, side))
            reinit!(fv, facet)
            for q_point in 1:getnquadpoints(fv)
                len += getdetJdV(fv, q_point)
                x = spatial_coordinate(fv, q_point, getcoordinates(facet))
                n = getnormal(fv, q_point)
                # On the circular arcs the normal must be (anti-)radial
                r = norm(x)
                if isapprox(r, r_in; atol = 1.0e-6) || isapprox(r, r_out; atol = 1.0e-6)
                    @test abs(n ⋅ (x / r)) ≈ 1.0 atol = 1.0e-6
                end
            end
        end
        @test len ≈ expected_lengths[side] rtol = 1.0e-10
    end
end

@testitem "FacetValues with a vector-valued field" begin
    using Ferrite

    geometry = createBSplineSquare(1.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 2)
    grid = IGAGrid{2}(geometry)

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    fqr = FacetQuadratureRule{RefQuadrilateral}(3)
    fv = FacetValues(fqr, ip^2, ip^2)

    dh = IGADofHandler(grid)
    add!(dh, :u, ip^2)
    close!(dh)

    # Assembling a constant traction must produce a total force = traction * edge length
    traction = Vec{2}((3.0, -2.0))
    f = zeros(ndofs(dh))
    fe = zeros(getnbasefunctions(fv))
    for facet in FacetIterator(dh, getfacetset(grid, "right"))
        reinit!(fv, facet)
        fill!(fe, 0.0)
        for q_point in 1:getnquadpoints(fv)
            dΓ = getdetJdV(fv, q_point)
            for i in 1:getnbasefunctions(fv)
                fe[i] += (traction ⋅ shape_value(fv, q_point, i)) * dΓ
            end
        end
        assemble!(f, celldofs(facet), fe)
    end
    @test sum(f[1:2:end]) ≈ traction[1]
    @test sum(f[2:2:end]) ≈ traction[2]
end
