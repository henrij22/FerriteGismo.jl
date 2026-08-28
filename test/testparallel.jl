@testitem "IGAInterpolation carries no per-cell state" begin
    geometry = createBSplineSquare(2.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 2)
    grid = IGAGrid{2}(geometry)
    @test getncells(grid) > 1

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))

    # No "current element" field to mutate or race on
    @test !ismutabletype(typeof(ip))
    @test fieldnames(typeof(ip)) == (:basis, :nbasefuns)

    # Ferrite's generic fallback `Base.copy(::Interpolation) = ip` is correct as-is now
    @test copy(ip) === ip

    ipVec = ip^2
    @test !ismutabletype(typeof(ipVec))
    @test copy(ipVec) === ipVec
end

@testitem "the same ip is safe to reinit! from independent CellValues at once" begin
    geometry = createBSplineSquare(2.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 2)
    grid = IGAGrid{2}(geometry)
    @test getncells(grid) > 1

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    qr = QuadratureRule{RefQuadrilateral}(2)

    dh = IGADofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)

    # Two independent CellValues built from the exact same ip, no copy anywhere
    cvA = CellValues(qr, ip, ip)
    cvB = CellValues(qr, ip, ip)
    @test Ferrite.function_interpolation(cvA) === Ferrite.function_interpolation(cvB) === ip

    ccA, ccB = CellCache(dh), CellCache(dh)
    reinit!(ccA, 1)
    reinit!(ccB, getncells(grid))
    reinit!(cvA, ccA)

    valsA_before = [shape_value(cvA, 1, i) for i in 1:getnbasefunctions(cvA)]
    detJdVA_before = getdetJdV(cvA, 1)

    # reinit!-ing cvB (sharing the same ip) to a different cell must not perturb cvA
    reinit!(cvB, ccB)

    valsA_after = [shape_value(cvA, 1, i) for i in 1:getnbasefunctions(cvA)]
    @test valsA_after == valsA_before
    @test getdetJdV(cvA, 1) == detJdVA_before

    # Sanity: cvB did pick up the other cell's (generally different) shape values
    valsB = [shape_value(cvB, 1, i) for i in 1:getnbasefunctions(cvB)]
    @test valsB != valsA_after || getdetJdV(cvB, 1) != detJdVA_before
end

@testitem "parallel assembly reproduces the serial stiffness matrix" begin
    using Base.Threads: nthreads, @threads
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

    # Parallel assembly: a pool of CellCache/CellValues (all sharing the same ip, no copy
    # needed) handed out via a Channel, smaller than getncells so scratch objects are
    # reused across tasks over the course of the loop.
    K_parallel = allocate_matrix(dh)
    assembler_lock = ReentrantLock()
    assembler_parallel = start_assemble(K_parallel)

    npool = 2 * nthreads() + 2
    pool = Channel{Tuple{CellCache, typeof(cv_serial), Matrix{Float64}}}(npool)
    for _ in 1:npool
        put!(pool, (CellCache(dh), CellValues(qr, ip, ip), zeros(n, n)))
    end

    @threads for cellid in 1:getncells(grid)
        cc, cv, ke_local = take!(pool)
        reinit!(cc, cellid)
        reinit!(cv, cc)
        assemble_element!(ke_local, cv)
        lock(assembler_lock) do
            return assemble!(assembler_parallel, celldofs(cc), ke_local)
        end
        put!(pool, (cc, cv, ke_local))
    end

    @test norm(K_parallel - K_serial) / norm(K_serial) < 1.0e-12
end
