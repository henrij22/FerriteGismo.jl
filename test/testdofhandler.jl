@testitem "DofHandler Test" begin
    geometry = createBSplineSquare(2.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 1)
    grid = IGAGrid{2}(geometry)
    basis = TinyGismo.basis(geometry)

    ip = IGAInterpolation{RefQuadrilateral}(basis)
    Ferrite.get_base_interpolation(ip) == ip

    dh = IGADofHandler(grid)
    add!(dh, :w, ip)
    close!(dh)

    @test ndofs_per_cell(dh) == getnbasefunctions(ip)
    @test ndofs(dh) == TinyGismo.size(basis)
    @test length(dh.cell_dofs) == 4 * getnbasefunctions(ip)
    @test length(dh.cell_dofs_offset) == 4
    @test dh.cell_dofs_offset[2] - dh.cell_dofs_offset[1] == getnbasefunctions(ip)
    @test dh.cell_dofs_offset[1] == 1
    @test celldofs(dh, 1) == dh.cell_dofs[dh.cell_dofs_offset[1]:(dh.cell_dofs_offset[2] - 1)]
    @test maximum(dh.cell_dofs) == TinyGismo.size(basis)

    for i in 1:4
        @test allunique(celldofs(dh, i))
    end

    @test length(dh.subdofhandlers) == 1
    @test first(dh.subdofhandlers).cellset == Ferrite.OrderedSet(1:4)

    dh = IGADofHandler(grid)
    add!(dh, :w, ip)
    add!(dh, :φ, ip)
    close!(dh)

    @test ndofs_per_cell(dh) == getnbasefunctions(ip) * 2
    @test ndofs(dh) == TinyGismo.size(basis) * 2
    @test length(dh.cell_dofs) == 4 * getnbasefunctions(ip) * 2
    @test length(dh.cell_dofs_offset) == 4
    @test dh.cell_dofs_offset[2] - dh.cell_dofs_offset[1] == getnbasefunctions(ip) * 2
    @test celldofs(dh, 1) == dh.cell_dofs[dh.cell_dofs_offset[1]:(dh.cell_dofs_offset[2] - 1)]
    @test maximum(dh.cell_dofs) == TinyGismo.size(basis) * 2

    for i in 1:4
        @test allunique(celldofs(dh, i))
    end

    ipVec = ip^1

    dhVec = IGADofHandler(grid)
    add!(dhVec, :w, ipVec)
    close!(dhVec)

    @test ndofs_per_cell(dhVec) == getnbasefunctions(ip)
    @test ndofs(dhVec) == TinyGismo.size(basis)
    @test length(dhVec.cell_dofs) == 4 * getnbasefunctions(ip)
    @test length(dhVec.cell_dofs_offset) == 4
    @test dhVec.cell_dofs_offset[2] - dhVec.cell_dofs_offset[1] == getnbasefunctions(ip)
    @test dhVec.cell_dofs_offset[1] == 1
    @test celldofs(dhVec, 1) == dhVec.cell_dofs[dhVec.cell_dofs_offset[1]:(dhVec.cell_dofs_offset[2] - 1)]
    @test maximum(dhVec.cell_dofs) == TinyGismo.size(basis)

    for i in 1:4
        @test allunique(celldofs(dhVec, i))
    end

    ipVec = ip^2
    @test Ferrite.n_components(ipVec) == 2
    Ferrite.get_base_interpolation(ipVec) == ip

    dhVec = IGADofHandler(grid)
    add!(dhVec, :w, ipVec)
    close!(dhVec)

    @test ndofs_per_cell(dhVec) == getnbasefunctions(ip) * 2
    @test ndofs(dhVec) == TinyGismo.size(basis) * 2
    @test dhVec.cell_dofs_offset[2] - dhVec.cell_dofs_offset[1] == getnbasefunctions(ip) * 2
    @test maximum(dhVec.cell_dofs) == TinyGismo.size(basis) * 2

    for i in 1:4
        @test allunique(celldofs(dhVec, i))
    end

    ipVec = ip^3

    dhVec = IGADofHandler(grid)
    add!(dhVec, :w, ipVec)
    close!(dhVec)

    @test ndofs_per_cell(dhVec) == getnbasefunctions(ip) * 3
    @test ndofs(dhVec) == TinyGismo.size(basis) * 3
    @test dhVec.cell_dofs_offset[2] - dhVec.cell_dofs_offset[1] == getnbasefunctions(ip) * 3
    @test maximum(dhVec.cell_dofs) == TinyGismo.size(basis) * 3

    for i in 1:4
        @test allunique(celldofs(dhVec, i))
    end
end
