@testitem "Vector-valued IGA interpolation" begin
    using Ferrite

    geometry = createBSplineSquare(1.0)
    degreeElevate!(geometry, 1)
    uniformRefine!(geometry, 2)
    grid = IGAGrid{2}(geometry)

    ip = IGAInterpolation{RefQuadrilateral}(TinyGismo.basis(geometry))
    qr = QuadratureRule{RefQuadrilateral}(2)

    cv_scalar = CellValues(qr, ip, ip^2)
    cv_vector = CellValues(qr, ip^2, ip^2)

    n_scalar = getnbasefunctions(cv_scalar)
    # A d-dimensional vector field has d times as many base functions.
    @test getnbasefunctions(cv_vector) == 2 * n_scalar

    e1 = Vec{2}((1.0, 0.0))
    e2 = Vec{2}((0.0, 1.0))

    for cell in CellIterator(grid)
        reinit!(cv_scalar, cell)
        reinit!(cv_vector, cell)

        for qp in 1:getnquadpoints(cv_vector)
            # The vectorized base functions are the scalar ones placed in a single
            # component (node-major, component-minor: u1x, u1y, u2x, u2y, ...).
            for i in 1:n_scalar
                φ = shape_value(cv_scalar, qp, i)
                ∇φ = shape_gradient(cv_scalar, qp, i)

                @test shape_value(cv_vector, qp, 2i - 1) ≈ φ * e1
                @test shape_value(cv_vector, qp, 2i) ≈ φ * e2

                # Gradient of a vector shape function is eᶜ ⊗ ∇φ.
                @test shape_gradient(cv_vector, qp, 2i - 1) ≈ e1 ⊗ ∇φ
                @test shape_gradient(cv_vector, qp, 2i) ≈ e2 ⊗ ∇φ
            end

            # Partition of unity: the values sum to (1, 1) per component ...
            val_sum = sum(shape_value(cv_vector, qp, i) for i in 1:getnbasefunctions(cv_vector))
            @test val_sum ≈ Vec{2}((1.0, 1.0))

            # ... and the gradients sum to zero.
            grad_sum = sum(shape_gradient(cv_vector, qp, i) for i in 1:getnbasefunctions(cv_vector))
            @test isapprox(grad_sum, zero(grad_sum); atol = 1.0e-10)
        end
    end
end
