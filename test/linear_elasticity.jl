# IGA adaptation of the Ferrite.jl linear elasticity tutorial:
# https://ferrite-fem.github.io/Ferrite.jl/stable/tutorials/linear_elasticity/
#
# We solve small-strain linear elasticity on a square, clamped on the left edge
# (cantilever) and loaded by a constant body force. The linear-elastic core
# (stiffness tensor `C` and the `∇N ⊡ C ⊡ ∇ˢʸᵐN` element stiffness) is identical to the
# Ferrite tutorial; the Neumann traction of the original example is replaced by a body
# force here, since facet-based integration is not yet available for IGA grids.

@testitem "Linear Elasticity Test" begin
    using Ferrite, SparseArrays

    # --- Geometry, grid, interpolation ---
    geometry = createBSplineSquare(1.0)
    degreeElevate!(geometry, 1)   # quadratic basis
    uniformRefine!(geometry, 3)
    grid = IGAGrid{2}(geometry)

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    qr = QuadratureRule{RefQuadrilateral}(3)
    cellvalues = CellValues(qr, ip^2, ip^2)

    dh = IGADofHandler(grid)
    add!(dh, :u, ip^2)
    close!(dh)

    # --- Boundary conditions: clamp the whole left edge (cantilever) ---
    ch = ConstraintHandler(dh)
    fixEdge!(dh, ch, :left, :u)
    close!(ch)

    # --- Linear elastic stiffness tensor (plane strain) ---
    Emod = 200.0e3
    ν = 0.3
    Gmod = Emod / (2 * (1 + ν))
    Kmod = Emod / (3 * (1 - 2ν))
    C = gradient(ϵ -> 2 * Gmod * dev(ϵ) + 3 * Kmod * vol(ϵ), zero(SymmetricTensor{2, 2}))

    # --- Constant body force (downward) ---
    b = Vec{2}((0.0, -1.0e3))

    function assemble_cell!(ke, fe, cellvalues, C, b)
        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)
            for i in 1:getnbasefunctions(cellvalues)
                Nᵢ = shape_value(cellvalues, q_point, i)
                ∇Nᵢ = shape_gradient(cellvalues, q_point, i)
                fe[i] += (b ⋅ Nᵢ) * dΩ
                for j in 1:getnbasefunctions(cellvalues)
                    ∇ˢʸᵐNⱼ = shape_symmetric_gradient(cellvalues, q_point, j)
                    ke[i, j] += (∇Nᵢ ⊡ C ⊡ ∇ˢʸᵐNⱼ) * dΩ
                end
            end
        end
        return ke, fe
    end

    function assemble_global(cellvalues, K, dh, C, b)
        n_basefuncs = getnbasefunctions(cellvalues)
        ke = zeros(n_basefuncs, n_basefuncs)
        fe = zeros(n_basefuncs)
        f = zeros(ndofs(dh))
        assembler = start_assemble(K, f)
        for cell in CellIterator(dh)
            reinit!(cellvalues, cell)
            fill!(ke, 0.0)
            fill!(fe, 0.0)
            assemble_cell!(ke, fe, cellvalues, C, b)
            assemble!(assembler, celldofs(cell), ke, fe)
        end
        return K, f
    end

    K = allocate_matrix(dh)
    K, f = assemble_global(cellvalues, K, dh, C, b)

    apply!(K, f, ch)
    u = K \ f

    # The clamped edge must remain fixed, and the free end should deflect downwards.
    uy = u[2:2:end]
    @test minimum(uy) < 0.0
    # Regression guard on the maximum (downward) tip deflection.
    @test minimum(uy) ≈ -0.014112710432389772 rtol = 1.0e-8
    @test maximum(abs.(u)) ≈ 0.014112710432389772 rtol = 1.0e-8

    mktempdir() do outdir
        VTKGridFile(joinpath(outdir, "linear_elasticity"), dh) do vtk
            write_solution(vtk, dh, u)
        end
    end
end
