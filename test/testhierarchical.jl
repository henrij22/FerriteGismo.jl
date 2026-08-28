# Hierarchical (HB/THB) patches.
#
# The defining difference from a tensor-product patch is that the number of basis functions
# acting on an element is not constant. FerriteGismo handles that by partitioning the cells
# into groups that share an active count and giving each group its own SubDofHandler, so most
# of what is checked here is that the grouping is a genuine partition and that the resulting
# global operator is still correct.

@testsnippet HierarchicalPatch begin
    # A locally refined unit square: quadratic, 4x4 coarse elements, with the lower-left
    # quadrant taken down two extra levels.
    function refinedSquare(; levels = 2, side = 1.0)
        geo = createBSplineSquare(side)
        degreeElevate!(geo)
        uniformRefine!(geo, 2)
        thb = THBSpline{2}(geo)
        for l in 2:levels
            refineElements!(thb, RefinementBox(l, 1:2, 1:2))
        end
        return thb
    end

    tensorSquare(; side = 1.0) = (geo = createBSplineSquare(side); degreeElevate!(geo); uniformRefine!(geo, 2); geo)

    function closedHandler(grid)
        dh = IGADofHandler(grid)
        for (cells, ip) in hierarchicalSubdomains(grid)
            sdh = SubDofHandler(dh, cells)
            add!(sdh, :u, ip)
        end
        close!(dh)
        return dh
    end

    # Assembles the mass and stiffness matrices of the Laplacian over the whole patch,
    # looping the subdomains as an application would.
    function assembleMassAndStiffness(dh; order = 4)
        M = allocate_matrix(dh)
        K = allocate_matrix(dh)
        aM, aK = start_assemble(M), start_assemble(K)
        qr = QuadratureRule{RefQuadrilateral}(order)
        area = 0.0
        for sdh in dh.subdofhandlers
            ip = only(sdh.field_interpolations)
            cv = CellValues(qr, ip, ip^2)
            n = getnbasefunctions(ip)
            Me, Ke = zeros(n, n), zeros(n, n)
            for cell in CellIterator(sdh)
                reinit!(cv, cell)
                fill!(Me, 0)
                fill!(Ke, 0)
                for q in 1:getnquadpoints(cv)
                    dΩ = getdetJdV(cv, q)
                    area += dΩ
                    for i in 1:n, j in 1:n
                        Me[i, j] += shape_value(cv, q, i) * shape_value(cv, q, j) * dΩ
                        Ke[i, j] += (shape_gradient(cv, q, i) ⋅ shape_gradient(cv, q, j)) * dΩ
                    end
                end
                assemble!(aM, celldofs(cell), Me)
                assemble!(aK, celldofs(cell), Ke)
            end
        end
        return M, K, area
    end
end

@testitem "Hierarchical: elements come from elementBoxes, not knotSpans" setup = [HierarchicalPatch] begin
    thb = refinedSquare()
    basis = TinyGismo.basis(thb)

    # TinyGismo refuses knotSpans for hierarchical bases, so the grid has to go the other way.
    @test_throws Exception knotSpans(basis)

    spans = FerriteGismo._elementSpans(basis, Val(2))
    @test length(spans) == Int(TinyGismo.numElements(basis))
    @test all(s -> all(s.lower .< s.upper), spans)
    @test all(s -> s.center ≈ (s.lower + s.upper) / 2, spans)
    # The elements partition the unit square
    @test sum(prod(s.upper - s.lower) for s in spans) ≈ 1.0

    # The tensor path is unchanged
    tensorSpans = FerriteGismo._elementSpans(TinyGismo.basis(tensorSquare()), Val(2))
    @test length(tensorSpans) == 9   # uniformRefine!(geo, 2) inserts two knots per direction
end

@testitem "Hierarchical: grid cells carry a varying number of control points" setup = [HierarchicalPatch] begin
    grid = IGAGrid{2}(refinedSquare())

    @test getncells(grid) == 12
    counts = [length(c.nodes) for c in grid.cells]
    # This is the whole reason subdomains are needed: the count is not constant.
    @test length(unique(counts)) > 1
    @test all(>=(9), counts)   # at least the (p+1)^2 of the coarse level
    @test all(1 .<= reduce(vcat, [c.nodes for c in grid.cells]) .<= length(grid.nodes))
end

@testitem "hierarchicalSubdomains partitions the cells by active count" setup = [HierarchicalPatch] begin
    grid = IGAGrid{2}(refinedSquare())
    groups = hierarchicalSubdomains(grid)

    @test length(groups) > 1
    cellsets = first.(groups)
    ips = last.(groups)

    # A genuine partition: every cell exactly once
    @test sort(reduce(vcat, cellsets)) == collect(1:getncells(grid))
    # Ascending active count, and each interpolation matches its group
    @test issorted(getnbasefunctions.(ips))
    for (cells, ip) in groups
        @test all(ci -> length(grid.cells[ci].nodes) == getnbasefunctions(ip), cells)
        @test Ferrite.getorder(ip) == 2
    end

    # A tensor grid has a single group, reproducing the plain interpolation
    tensorGrid = IGAGrid{2}(tensorSquare())
    tensorGroups = hierarchicalSubdomains(tensorGrid)
    @test length(tensorGroups) == 1
    @test first(tensorGroups)[1] == collect(1:getncells(tensorGrid))
    @test getnbasefunctions(first(tensorGroups)[2]) ==
        getnbasefunctions(IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(tensorSquare())))
end

@testitem "IGAInterpolation rejects a hierarchical basis and points at the fix" setup = [HierarchicalPatch] begin
    basis = TinyGismo.basis(refinedSquare())
    err = try
        IGAInterpolation{RefQuadrilateral}(basis)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("hierarchicalSubdomains", err.msg)

    # The explicit form is what the grouping uses
    ip = IGAInterpolation{RefQuadrilateral}(basis, 12)
    @test getnbasefunctions(ip) == 12
    @test Ferrite.getorder(ip) == 2
end

@testitem "Hierarchical: dofs are numbered globally over the control points" setup = [HierarchicalPatch] begin
    thb = refinedSquare()
    grid = IGAGrid{2}(thb)
    dh = closedHandler(grid)

    # One dof per control point of the whole patch, shared across subdomains
    @test ndofs(dh) == Int(size(TinyGismo.basis(thb)))
    @test length(dh.subdofhandlers) == length(hierarchicalSubdomains(grid))

    # Per-cell dof counts follow the subdomain, and celldofs agrees with the cell's actives
    for cell in CellIterator(dh)
        ci = cellid(cell)
        @test Ferrite.ndofs_per_cell(dh, ci) == length(grid.cells[ci].nodes)
        @test length(celldofs(cell)) == length(grid.cells[ci].nodes)
        @test sort(celldofs(cell)) == sort(grid.cells[ci].nodes)
    end

    # Every dof is reached by at least one cell
    @test sort(unique(dh.cell_dofs)) == collect(1:ndofs(dh))
end

@testitem "Hierarchical: the assembled operator is exact" setup = [HierarchicalPatch] begin
    using LinearAlgebra
    thb = refinedSquare(; levels = 3, side = 2.0)   # [0,2]^2, area 4
    dh = closedHandler(IGAGrid{2}(thb))
    M, K, area = assembleMassAndStiffness(dh)

    @test area ≈ 4.0
    # sum(M) = ∫(ΣNᵢ)(ΣNⱼ) = ∫1 = area. This holds only if the shape values are a partition
    # of unity *and* the per-cell dofs scatter into the right global entries, so it checks
    # the subdomain bookkeeping as much as the basis.
    @test sum(M) ≈ 4.0
    # Constants are in the kernel of the Laplacian.
    @test norm(K * ones(ndofs(dh))) < 1.0e-10
    @test issymmetric(Matrix(K))
end

@testitem "Hierarchical: an unrefined THB patch reproduces the tensor patch" setup = [HierarchicalPatch] begin
    using LinearAlgebra
    # With a single level a THB basis *is* the tensor basis, so the whole pipeline must
    # produce the identical system -- one group, same dofs, same matrices.
    tensorGeo = tensorSquare()
    thb = THBSpline{2}(tensorGeo)

    thbGrid, tensorGrid = IGAGrid{2}(thb), IGAGrid{2}(tensorGeo)
    @test getncells(thbGrid) == getncells(tensorGrid)

    thbDh = closedHandler(thbGrid)
    @test length(thbDh.subdofhandlers) == 1

    tensorDh = IGADofHandler(tensorGrid)
    add!(tensorDh, :u, IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(tensorGeo)))
    close!(tensorDh)

    @test ndofs(thbDh) == ndofs(tensorDh)
    @test thbDh.cell_dofs == tensorDh.cell_dofs
    @test thbDh.cell_dofs_offset == tensorDh.cell_dofs_offset

    Mt, Kt, _ = assembleMassAndStiffness(thbDh)
    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(tensorGeo))
    Mr, Kr = allocate_matrix(tensorDh), allocate_matrix(tensorDh)
    aM, aK = start_assemble(Mr), start_assemble(Kr)
    qr = QuadratureRule{RefQuadrilateral}(4)
    cv = CellValues(qr, ip, ip^2)
    n = getnbasefunctions(ip)
    Me, Ke = zeros(n, n), zeros(n, n)
    for cell in CellIterator(tensorDh)
        reinit!(cv, cell)
        fill!(Me, 0)
        fill!(Ke, 0)
        for q in 1:getnquadpoints(cv)
            dΩ = getdetJdV(cv, q)
            for i in 1:n, j in 1:n
                Me[i, j] += shape_value(cv, q, i) * shape_value(cv, q, j) * dΩ
                Ke[i, j] += (shape_gradient(cv, q, i) ⋅ shape_gradient(cv, q, j)) * dΩ
            end
        end
        assemble!(aM, celldofs(cell), Me)
        assemble!(aK, celldofs(cell), Ke)
    end
    @test Matrix(Mt) ≈ Matrix(Mr)
    @test Matrix(Kt) ≈ Matrix(Kr)
end

@testitem "Hierarchical: boundary conditions work across subdomains" setup = [HierarchicalPatch] begin
    grid = IGAGrid{2}(refinedSquare())
    dh = closedHandler(grid)

    ch = ConstraintHandler(dh)
    for side in (:left, :right, :bottom, :top)
        fixEdge!(dh, ch, side, :u)
    end
    close!(ch)

    @test !isempty(ch.prescribed_dofs)
    @test all(1 .<= ch.prescribed_dofs .<= ndofs(dh))
    @test length(unique(ch.prescribed_dofs)) == length(ch.prescribed_dofs)

    # edgeControlPoints resolves through any subdomain, not just a single one
    @test !isempty(edgeControlPoints(dh, :left, :u))
end

@testitem "Hierarchical: Poisson solve converges to the analytic value" setup = [HierarchicalPatch] begin
    # -Δu = 1 on the unit square with u = 0 on the boundary has u(1/2,1/2) = 0.0736713...
    function solveCentre(levels)
        grid = IGAGrid{2}(refinedSquare(; levels))
        dh = closedHandler(grid)
        ch = ConstraintHandler(dh)
        for side in (:left, :right, :bottom, :top)
            fixEdge!(dh, ch, side, :u)
        end
        close!(ch)

        K = allocate_matrix(dh)
        f = zeros(ndofs(dh))
        assembler = start_assemble(K, f)
        qr = QuadratureRule{RefQuadrilateral}(4)
        for sdh in dh.subdofhandlers
            ip = only(sdh.field_interpolations)
            cv = CellValues(qr, ip, ip^2)
            n = getnbasefunctions(ip)
            Ke, fe = zeros(n, n), zeros(n)
            for cell in CellIterator(sdh)
                reinit!(cv, cell)
                fill!(Ke, 0)
                fill!(fe, 0)
                for q in 1:getnquadpoints(cv)
                    dΩ = getdetJdV(cv, q)
                    for i in 1:n
                        fe[i] += shape_value(cv, q, i) * dΩ
                        for j in 1:n
                            Ke[i, j] += (shape_gradient(cv, q, i) ⋅ shape_gradient(cv, q, j)) * dΩ
                        end
                    end
                end
                assemble!(assembler, celldofs(cell), Ke, fe)
            end
        end
        apply!(K, f, ch)
        u = K \ f
        return FerriteGismo.interpolate(FerriteGismo._fieldInterpolation(dh, :u), u, [0.5, 0.5])
    end

    exact = 0.07367135328
    e2, e3 = abs(solveCentre(2) - exact), abs(solveCentre(3) - exact)
    @test e2 < 1.0e-3
    @test e3 < e2          # refining further gets closer
end

@testitem "Hierarchical: export evaluates the field correctly" setup = [HierarchicalPatch] begin
    grid = IGAGrid{2}(refinedSquare())
    dh = closedHandler(grid)

    # A coefficient vector of ones is the constant one function, by partition of unity.
    vals = evaluateAtExportNodes(dh, ones(ndofs(dh)), :u)
    @test all(v -> isapprox(v, 1.0; atol = 1.0e-12), vals)

    # The export mesh is built element by element, so it draws the hierarchical layout
    # exactly rather than a uniform lattice over it.
    @test isHierarchical(grid)
    @test getncells(exportGrid(grid)) == getncells(grid)
    @test getncells(parameterSpaceGrid(grid)) == getncells(grid)
    @test getncells(parameterSpaceGrid(grid; subdivision = 3)) == 9 * getncells(grid)

    # Its parametric nodes are the element corners, and they cover the unit square.
    ps = parameterSpaceGrid(grid)
    coords = [n.x for n in ps.nodes]
    @test length(coords) == 4 * getncells(grid)      # corners are duplicated per element
    @test all(c -> all(0.0 .<= c .<= 1.0), coords)

    # Notions that presuppose a tensor lattice are refused rather than answered wrongly.
    @test_throws Exception breakpoints(grid, 1)
    @test_throws Exception numElementsPerDirection(grid, 1)

    # The tensor path keeps its lattice mesh, sharing nodes between cells.
    tensorGrid = IGAGrid{2}(tensorSquare())
    @test !isHierarchical(tensorGrid)
    @test getncells(parameterSpaceGrid(tensorGrid)) == getncells(tensorGrid)
    @test length(parameterSpaceGrid(tensorGrid).nodes) == 16   # 4x4 lattice of a 3x3 mesh
end
