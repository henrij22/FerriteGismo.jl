# See https://ferrite-fem.github.io/Ferrite.jl/stable/tutorials/heat_equation/
# for explanations of this tutorial

@testitem "Heat Equation Test" begin
    using Ferrite, SparseArrays

    geometry = createBSplineSquare(2.0)
    uniformRefine!(geometry, 19)
    grid = IGAGrid{2}(geometry)

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    qr = QuadratureRule{RefQuadrilateral}(2)
    cellvalues = CellValues(qr, ip, ip^2)

    dh = IGADofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)

    K = allocate_matrix(dh)

    ch = ConstraintHandler(dh)
    fixEdge!(dh, ch, :left, :u)
    fixEdge!(dh, ch, :right, :u)
    fixEdge!(dh, ch, :bottom, :u)
    fixEdge!(dh, ch, :top, :u)

    close!(ch)

    function assemble_element!(Ke::Matrix, fe::Vector, cellvalues::CellValues)
        n_basefuncs = getnbasefunctions(cellvalues)
        # Reset to 0
        fill!(Ke, 0)
        fill!(fe, 0)
        # Loop over quadrature points
        for q_point in 1:getnquadpoints(cellvalues)
            # Get the quadrature weight
            dΩ = getdetJdV(cellvalues, q_point)
            # Loop over test shape functions
            for i in 1:n_basefuncs
                δu = shape_value(cellvalues, q_point, i)
                ∇δu = shape_gradient(cellvalues, q_point, i)
                # Add contribution to fe
                fe[i] += δu * dΩ
                # Loop over trial shape functions
                for j in 1:n_basefuncs
                    ∇u = shape_gradient(cellvalues, q_point, j)
                    # Add contribution to Ke
                    Ke[i, j] += (∇δu ⋅ ∇u) * dΩ
                end
            end
        end
        return Ke, fe
    end

    function assemble_global(cellvalues::CellValues, K::SparseMatrixCSC, dh::IGADofHandler)
        # Allocate the element stiffness matrix and element force vector
        n_basefuncs = getnbasefunctions(cellvalues)
        Ke = zeros(n_basefuncs, n_basefuncs)
        fe = zeros(n_basefuncs)
        # Allocate global force vector f
        f = zeros(ndofs(dh))
        # Create an assembler
        assembler = start_assemble(K, f)
        # Loop over all cells
        for cell in CellIterator(dh)
            # Reinitialize cellvalues for this cell
            reinit!(cellvalues, cell)
            # Compute element contribution
            assemble_element!(Ke, fe, cellvalues)
            # Assemble Ke and fe into K and f
            assemble!(assembler, celldofs(cell), Ke, fe)
        end
        return K, f
    end

    K, f = assemble_global(cellvalues, K, dh)

    apply!(K, f, ch)
    u = K \ f

    @test maximum(u) ≈ 0.295267863770736

    VTKGridFile("heat_equation", dh) do vtk
        write_solution(vtk, dh, u)
    end
end
