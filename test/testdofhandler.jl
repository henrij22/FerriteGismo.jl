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

@testitem "Several SubDofHandlers share one dof numbering per field" begin
    # Subdomains partition the *cells* of one patch, so a field is numbered once over the
    # control points and every subdomain that carries it sees the same dofs.
    geometry = createBSplineSquare(1.0)
    uniformRefine!(geometry, 2)
    grid = IGAGrid{2}(geometry)
    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))

    n = getncells(grid)
    lower, upper = collect(1:(n ÷ 2)), collect((n ÷ 2 + 1):n)

    split = IGADofHandler(grid)
    for cells in (lower, upper)
        sdh = SubDofHandler(split, cells)
        add!(sdh, :u, ip)
        add!(sdh, :v, ip)
    end
    close!(split)

    whole = IGADofHandler(grid)
    add!(whole, :u, ip)
    add!(whole, :v, ip)
    close!(whole)

    # Splitting the cells changes nothing about the dofs
    @test ndofs(split) == ndofs(whole)
    @test split.cell_dofs == whole.cell_dofs
    @test split.cell_dofs_offset == whole.cell_dofs_offset
    @test fieldOffset(split, :u) == fieldOffset(whole, :u) == 0
    @test fieldOffset(split, :v) == fieldOffset(whole, :v) == size(ip.basis)

    # Every cell is claimed by exactly one subdomain
    @test sort(vcat(collect.(getproperty.(split.subdofhandlers, :cellset))...)) == collect(1:n)
end

@testitem "fieldOffset is the global offset, not a SubDofHandler-local one" begin
    geometry = createBSplineSquare(1.0)
    uniformRefine!(geometry, 2)
    grid = IGAGrid{2}(geometry)
    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    n = getncells(grid)

    # :a lives only on the first subdomain, :b only on the second, so :b is field 1 of its
    # own SubDofHandler while being the second field globally.
    dh = IGADofHandler(grid)
    sdhA = SubDofHandler(dh, collect(1:(n ÷ 2)))
    add!(sdhA, :a, ip)
    sdhB = SubDofHandler(dh, collect((n ÷ 2 + 1):n))
    add!(sdhB, :b, ip)
    close!(dh)

    @test Ferrite._find_field(sdhB, :b) == 1        # local index
    @test fieldOffset(dh, :a) == 0
    @test fieldOffset(dh, :b) == size(ip.basis)     # not 0
    @test ndofs(dh) == 2 * size(ip.basis)
    @test_throws Exception fieldOffset(dh, :nope)
end

@testitem "A field must agree across the SubDofHandlers that carry it" begin
    coarse = createBSplineSquare(1.0)
    uniformRefine!(coarse, 2)
    grid = IGAGrid{2}(coarse)

    fine = createBSplineSquare(1.0)
    uniformRefine!(fine, 2)
    degreeElevate!(fine)

    ipCoarse = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(coarse))
    ipFine = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(fine))
    n = getncells(grid)

    # Two subdomains describing :u by different spline spaces would silently corrupt the
    # shared numbering, so closing has to reject it. (Ferrite already rejects a differing
    # *component* count at `add!` time; the basis behind an IGAInterpolation is opaque to it,
    # which is the gap this closes.)
    dh = IGADofHandler(grid)
    add!(SubDofHandler(dh, collect(1:(n ÷ 2))), :u, ipCoarse)
    add!(SubDofHandler(dh, collect((n ÷ 2 + 1):n)), :u, ipFine)
    @test size(ipCoarse.basis) != size(ipFine.basis)
    @test_throws "different interpolations" close!(dh)
end
