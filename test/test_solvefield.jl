@testset "solveField: backends and reconstruction" begin
    with_test_gmsh() do
        @testset "RHS handling and linear solvers" begin
            structured_box_mesh(n=2)

            mat = Material("body")
            U = Field(
                [mat],
                type=:VectorField,
                dim=3,
                fieldName=:u
            )

            K = ∫(ε(U) ⋅ D(:Solid, mat) ⋅ ε(U))
            f1 = ∫(U ⋅ [1.0, 0.0, 0.0], Γ="right")
            f2 = 2.0 * f1

            support = [
                BoundaryCondition("left", ux=0.0, uy=0.0, uz=0.0)
            ]

            # Multiple right-hand sides must reproduce separate solves.
            f_multi = VectorField(
                Matrix{Float64}[],
                hcat(f1.a[:, 1], f2.a[:, 1]),
                [0.0, 1.0],
                Int[],
                2,
                f1.type,
                f1.model
            )

            u_multi = solveField(K, f_multi; support=support)
            u1 = solveField(K, f1; support=support)
            u2 = solveField(K, f2; support=support)

            @test size(u_multi.a, 2) == 2
            @test relerr(u_multi.a[:, 1], u1.a[:, 1]) < 1e-12
            @test relerr(u_multi.a[:, 2], u2.a[:, 1]) < 1e-12

            # Direct solvers and symmetric wrapper.
            u_lu = solveField(
                Symmetric(K),
                f1;
                support=support,
                solver=:lu
            )

            u_chol = solveField(
                Symmetric(K),
                f1;
                support=support,
                solver=:cholesky
            )

            u_auto = solveField(
                Symmetric(K),
                f1;
                support=support
            )

            @test relerr(u_lu.a, u1.a) < 1e-10
            @test relerr(u_chol.a, u1.a) < 1e-10
            @test relerr(u_auto.a, u_chol.a) < 1e-12

            # Iterative solvers.
            u_cg = solveField(
                Symmetric(K),
                f1;
                support=support,
                solver=:cg,
                solveroptions=(reltol=1e-10, maxiter=size(K.A, 1))
            )

            u_gmres = solveField(
                K,
                f1;
                support=support,
                solver=:gmres,
                solveroptions=(reltol=1e-10,maxiter=size(K.A, 1))
            )

            @test relerr(u_cg.a, u1.a) < 1e-7
            @test relerr(u_gmres.a, u1.a) < 1e-7

        end

        @testset "Reduced-order interpolation" begin
            structured_box_mesh(n=1, order=2)

            mat = Material("body")
            U = Field(
                [mat],
                type=:VectorField,
                dim=3,
                fieldName=:u,
                reducedOrder=true
            )

            K = ∫(ε(U) ⋅ D(:Solid, mat) ⋅ ε(U))
            f = ∫(U ⋅ [1.0, 0.0, 0.0], Γ="right")

            support = [
                BoundaryCondition("left", ux=0.01, uy=0.0, uz=0.0)
            ]

            u = solveField(K, f; support=support)
            u_sym = solveField(Symmetric(K), f; support=support)

            T, R = reductionMatrices(U)

            @test size(T, 1) == size(K.A, 1)
            @test size(T, 2) < size(T, 1)
            @test size(R) == reverse(size(T))
            @test relerr(u_sym.a, u.a) < 1e-10
            @test u.a ≈ T * (R * u.a) atol=1e-10
            left_x = constrainedDoFs(U, [BoundaryCondition("left", ux=0.0)])
            @test u.a[left_x, 1] ≈ fill(0.01, length(left_x)) atol=1e-12

        end
    end
end

@testset "solveField: analytical static solutions" begin
    with_test_gmsh() do
        line_mesh(n=2, order=2)
        P = Field([Material("body")], type=:ScalarField, dim=1, fieldName=:p)
        K = ∫(Grad(P) ⋅ Grad(P))
        f = ∫(P ⋅ 2.0)
        support = [BoundaryCondition("left"; field=P, p=1.0),
                   BoundaryCondition("right"; field=P, p=2.0)]
        # -p''=2, p(0)=1, p(1)=2 gives p=1+2x-x², exactly in P2.
        exact = scalarField(P, "body", (x,y,z)->1+2x-x^2)
        u = solveField(K, f; support=support)
        @test u isa ScalarField
        @test u.a ≈ exact.a rtol=1e-11 atol=1e-12
        free, fixed = freeDoFs(P, support), constrainedDoFs(P, support)
        @test sort(vcat(free, fixed)) == collect(1:ndofs(P))
        @test isempty(intersect(free, fixed))
        @test norm((K.A * u.a - f.a)[free, :]) < 1e-10
    end

    with_test_gmsh() do
        structured_rect_mesh(n=2, order=1)
        mat = Material("body")
        U = Field([mat], type=:VectorField, dim=2, fieldName=:u)
        K = ∫(SymGrad(U) ⋅ D(:PlaneStress, mat) ⋅ SymGrad(U))
        f = ∫(U ⋅ [1.0, 0.0]; Γ="right")
        support = [BoundaryCondition("left"; field=U, ux=0.0),
                   BoundaryCondition("bottom"; field=U, uy=0.0)]
        u = solveField(Symmetric(K), f; support=support)
        exact = vectorField(U, "body", [(x,y,z)->x/mat.E,
                                       (x,y,z)->-mat.ν*y/mat.E, 0.0])
        @test relerr(u.a, exact.a) < 1e-10
        @test probe(u, 1.0, 1.0, 0.0)[1:2] ≈ [1/mat.E, -mat.ν/mat.E] rtol=1e-10
        free = freeDoFs(U, support)
        residual = K.A * u.a[:, 1] - f.a[:, 1]
        @test norm(residual[free]) / norm(f.a[free, 1]) < 1e-10
        @test all(iszero, u.a[constrainedDoFs(U, support), :])
    end
end


@testset "Shared linear backend: cached solves and sparse RHS" begin
    A=sparse([3.0 -1.0; -1.0 2.0])
    B=[1.0 2.0; 3.0 4.0]
    for solver in (:backslash,:lu,:lu_no_ordering,:cholesky,:qr,:cg,:gmres)
        options=solver in (:cg,:gmres) ? (reltol=1e-12,) : (;)
        prepared=LowLevelFEM.prepare_linear_solver(Symmetric(A);solver,solveroptions=options)
        @test LowLevelFEM.solve_prepared_linear_system(prepared,B) ≈ A\B
        for j in axes(B,2)
            b=sparsevec(B[:,j])
            @test LowLevelFEM.solve_prepared_linear_system(prepared,b) ≈ A\B[:,j]
        end
    end
    rectangular=sparse([1.0 0.0; 0.0 2.0; 1.0 1.0])
    rhs=reshape([1.0,2.0,3.0],:,1)
    @test LowLevelFEM.solve_linear_system(rectangular,rhs) ≈ rectangular\rhs
end
