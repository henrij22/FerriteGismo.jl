@testitem "copy(::IGAInterpolation) is independent of the original" begin
    geometry = createBSplineSquare(2.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 2)
    grid = IGAGrid{2}(geometry)

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    ipCopy = copy(ip)

    # A genuinely different object, so mutating one's `currentElement` (as `reinit!` does)
    # cannot leak into the other.
    @test ipCopy !== ip
    @test ipCopy.currentElement == ip.currentElement

    # The read-only pieces are still shared (no need to duplicate the G+Smo basis itself).
    @test ipCopy.basis === ip.basis
    @test ipCopy.nbasefuns == ip.nbasefuns

    ipCopy.currentElement = grid.knotSpans[end]
    @test ip.currentElement == grid.knotSpans[1]
    @test ipCopy.currentElement == grid.knotSpans[end]

    # `copy` must also reach through `VectorizedInterpolation` (`ip^2`), since vector-valued
    # fields are the common case.
    ipVec = ip^2
    ipVecCopy = copy(ipVec)
    @test ipVecCopy !== ipVec
    @test ipVecCopy.ip !== ipVec.ip
    @test ipVecCopy.ip.basis === ipVec.ip.basis

    ipVecCopy.ip.currentElement = grid.knotSpans[end]
    @test ipVec.ip.currentElement == grid.knotSpans[1]
end

@testitem "copy(::CellValues) gives task-local IGA interpolations" begin
    geometry = createBSplineSquare(2.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 2)
    grid = IGAGrid{2}(geometry)
    @test getncells(grid) > 1

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    qr = QuadratureRule{RefQuadrilateral}(2)
    cellvalues = CellValues(qr, ip, ip)

    dh = IGADofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)

    cvA = copy(cellvalues)
    cvB = copy(cellvalues)
    @test Ferrite.function_interpolation(cvA) !== Ferrite.function_interpolation(cvB)

    ccA, ccB = CellCache(dh), CellCache(dh)
    reinit!(ccA, 1)
    reinit!(ccB, getncells(grid))
    reinit!(cvA, ccA)

    # Snapshot cvA's shape values before touching cvB, so that a stale reference to a
    # shared `currentElement` (the bug this guards against) would corrupt them.
    valsA_before = [shape_value(cvA, 1, i) for i in 1:getnbasefunctions(cvA)]
    detJdVA_before = getdetJdV(cvA, 1)

    # reinit!-ing the *other* copy to a different cell must not perturb cvA at all.
    reinit!(cvB, ccB)

    valsA_after = [shape_value(cvA, 1, i) for i in 1:getnbasefunctions(cvA)]
    @test valsA_after == valsA_before
    @test getdetJdV(cvA, 1) == detJdVA_before

    # Sanity: cvB really did pick up the other cell's (generally different) shape values.
    valsB = [shape_value(cvB, 1, i) for i in 1:getnbasefunctions(cvB)]
    @test valsB != valsA_after || getdetJdV(cvB, 1) != detJdVA_before
end

@testitem "parallel assembly reproduces the serial stiffness matrix" begin
    using Base.Threads: nthreads, threadid, @threads
    using LinearAlgebra: norm

    geometry = createBSplineSquare(1.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 3)
    grid = IGAGrid{2}(geometry)

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    qr = QuadratureRule{RefQuadrilateral}(3)

    dh = IGADofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)

    n = ndofs_per_cell(dh)

    function assemble_element!(ke, cv)
        fill!(ke, 0.0)
        for q_point in 1:getnquadpoints(cv)
            dΩ = getdetJdV(cv, q_point)
            for i in 1:getnbasefunctions(cv)
                ∇φi = shape_gradient(cv, q_point, i)
                for j in 1:getnbasefunctions(cv)
                    ∇φj = shape_gradient(cv, q_point, j)
                    ke[i, j] += (∇φi ⋅ ∇φj) * dΩ
                end
            end
        end
        return ke
    end

    # Serial reference
    K_serial = allocate_matrix(dh)
    assembler_serial = start_assemble(K_serial)
    cv_serial = CellValues(qr, ip, ip)
    ke = zeros(n, n)
    for cell in CellIterator(dh)
        reinit!(cv_serial, cell)
        assemble_element!(ke, cv_serial)
        assemble!(assembler_serial, celldofs(cell), ke)
    end

    # "Parallel" assembly: every physical thread works off its own `copy` of `CellValues`/
    # `CellCache` (the pattern the docs recommend), indexed by `threadid()` under `:static`
    # scheduling so this is a real one-scratch-set-per-thread setup, not just an aliased
    # re-use of the same objects. This exercises real concurrency whenever the test process
    # has more than one thread.
    K_parallel = allocate_matrix(dh)
    assembler_lock = ReentrantLock()
    assembler_parallel = start_assemble(K_parallel)

    ccs = [CellCache(dh) for _ in 1:nthreads()]
    cvs = [copy(cv_serial) for _ in 1:nthreads()]
    kes = [zeros(n, n) for _ in 1:nthreads()]

    @threads :static for cellid in 1:getncells(grid)
        tid = threadid()
        cc, cv, ke_local = ccs[tid], cvs[tid], kes[tid]
        reinit!(cc, cellid)
        reinit!(cv, cc)
        assemble_element!(ke_local, cv)
        lock(assembler_lock) do
            return assemble!(assembler_parallel, celldofs(cc), ke_local)
        end
    end

    @test norm(K_parallel - K_serial) / norm(K_serial) < 1.0e-12
end
