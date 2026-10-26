@testset "Legacy stiffness and DSL equivalence" begin
    with_test_gmsh() do
        structured_box_mesh(n=2, order=1)
        mat = Material("body")
        Plegacy = Problem([mat])
        U = Field([mat], type=:VectorField, dim=3, fieldName=:u)
        Klegacy = stiffnessMatrix(Plegacy)
        Kdsl = ∫(SymGrad(U) ⋅ D(:Solid, mat) ⋅ SymGrad(U))
        @test size(Kdsl.A) == size(Klegacy.A)
        @test norm(Kdsl.A - Klegacy.A) / norm(Klegacy.A) < 1e-12

        f = ∫(U ⋅ [1.0, 0.0, 0.0]; Γ="right")
        support = [BoundaryCondition("left", ux=0.0, uy=0.0, uz=0.0)]
        u = solveField(Kdsl, f; support=support)
        u_cg = solveField(Symmetric(Kdsl), f; support=support, solver=:cg,
                         solveroptions=(reltol=1e-10, maxiter=size(Kdsl.A, 1)))
        u_old = solveField(Symmetric(Kdsl), f; support=support, iterative=true,
                          solveroptions=(reltol=1e-10, maxiter=size(Kdsl.A, 1)))
        u_order = solveField(Kdsl, f; support=support, ordering=false)
        @test relerr(u_old.a, u_cg.a) < 1e-12
        @test relerr(u_order.a, u.a) < 1e-10

        # Preserve the raw-result compatibility wrapper from the eigen dev tests.
        M = ∫(U ⋅ U)
        oldλ, oldϕ = solveEigenProblem(Kdsl, M; n=3, directSolver=true, support=support)
        λ, ϕ = solveEigen(Kdsl, M; n=3, solver=:direct, output=:raw, support=support)
        @test oldλ ≈ λ
        test_eigen_residual(Kdsl.A, M.A, oldλ, oldϕ;
                            transform=sparse(I, ndofs(U), ndofs(U))[:, freeDoFs(U, support)])
    end
end
