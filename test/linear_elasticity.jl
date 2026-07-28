# IGA adaptation of the Ferrite.jl linear elasticity tutorial:
# https://ferrite-fem.github.io/Ferrite.jl/stable/tutorials/linear_elasticity/
#
# We solve small-strain linear elasticity on a square, clamped on the left edge
# (cantilever) and loaded by a constant Neumann traction on the right edge. The
# linear-elastic core (stiffness tensor `C` and the `∇N ⊡ C ⊡ ∇ˢʸᵐN` element stiffness) is
# identical to the Ferrite tutorial; the traction is integrated with `FacetValues` over the
# boundary facet set of the right edge.

@testitem "Linear Elasticity Test" begin
    using Ferrite, SparseArrays

    # --- Geometry, grid, interpolation ---
    geometry = createBSplineSquare(1.0)
    degreeElevate!(geometry, 1)   # quadratic basis
    uniformRefine!(geometry, 3)
    grid = IGAGrid{2}(geometry)

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    qr = QuadratureRule{RefQuadrilateral}(3)
    fqr = FacetQuadratureRule{RefQuadrilateral}(3)
    cellvalues = CellValues(qr, ip^2, ip^2)
    facetvalues = FacetValues(fqr, ip^2, ip^2)

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

    function assemble_cell!(ke, cellvalues, C)
        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)
            for i in 1:getnbasefunctions(cellvalues)
                ∇Nᵢ = shape_gradient(cellvalues, q_point, i)
                for j in 1:getnbasefunctions(cellvalues)
                    ∇ˢʸᵐNⱼ = shape_symmetric_gradient(cellvalues, q_point, j)
                    ke[i, j] += (∇Nᵢ ⊡ C ⊡ ∇ˢʸᵐNⱼ) * dΩ
                end
            end
        end
        return ke
    end

    function assemble_global!(K, dh, cellvalues, C)
        n_basefuncs = getnbasefunctions(cellvalues)
        ke = zeros(n_basefuncs, n_basefuncs)
        assembler = start_assemble(K)
        for cell in CellIterator(dh)
            reinit!(cellvalues, cell)
            fill!(ke, 0.0)
            assemble_cell!(ke, cellvalues, C)
            assemble!(assembler, celldofs(cell), ke)
        end
        return K
    end

    # --- Neumann traction on the right edge ---
    function assemble_external_forces!(f_ext, dh, facetset, facetvalues, traction)
        fe_ext = zeros(getnbasefunctions(facetvalues))
        for facet in FacetIterator(dh, facetset)
            reinit!(facetvalues, facet)
            fill!(fe_ext, 0.0)
            for q_point in 1:getnquadpoints(facetvalues)
                dΓ = getdetJdV(facetvalues, q_point)
                for i in 1:getnbasefunctions(facetvalues)
                    Nᵢ = shape_value(facetvalues, q_point, i)
                    fe_ext[i] += (traction ⋅ Nᵢ) * dΓ
                end
            end
            assemble!(f_ext, celldofs(facet), fe_ext)
        end
        return f_ext
    end

    K = allocate_matrix(dh)
    assemble_global!(K, dh, cellvalues, C)

    traction = Vec{2}((0.0, -1.0e3))
    f = zeros(ndofs(dh))
    assemble_external_forces!(f, dh, getfacetset(grid, "right"), facetvalues, traction)

    # The total applied force must equal traction × edge length.
    @test sum(f[2:2:end]) ≈ -1.0e3
    @test sum(f[1:2:end]) ≈ 0.0 atol = 1.0e-12

    apply!(K, f, ch)
    u = K \ f

    # The clamped edge must remain fixed, and the free end should deflect downwards.
    uy = u[2:2:end]
    @test minimum(uy) < 0.0
    # Regression guard on the maximum (downward) tip deflection.
    @test minimum(uy) ≈ -0.03392555078206797 rtol = 1.0e-8
    @test maximum(abs.(u)) ≈ 0.03392555078206797 rtol = 1.0e-8

    mktempdir() do outdir
        VTKGridFile(joinpath(outdir, "linear_elasticity"), dh) do vtk
            write_solution(vtk, dh, u)
        end
    end
end
