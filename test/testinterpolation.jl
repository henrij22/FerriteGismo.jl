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
