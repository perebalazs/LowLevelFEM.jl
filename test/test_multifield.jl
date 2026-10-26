@testset "Multifield block systems" begin
    with_test_gmsh() do
        (; P, V, A, B, C, D, K, fv_x, fp, F_x, pres, supp_top, supp_bottom, support) = stokes_fixture()
        @testset "One-field block return type" begin
            Kblock = SystemMatrix(reshape([A + D], 1, 1))
            Fblock = SystemVector([fv_x])

            v_block = solveField(
                Kblock,
                Fblock;
                support=[supp_top, supp_bottom]
            )

            v_single = solveField(
                A + D,
                fv_x;
                support=[supp_top, supp_bottom]
            )

            @test v_block isa VectorField
            @test !(v_block isa Tuple)
            @test relerr(v_block.a, v_single.a) < 1e-12
        end

        @testset "Reduced multifield system and symmetric fallback" begin
            v, p = solveField(
                K,
                F_x;
                support=support
            )

            prepared = LowLevelFEM.prepare_multifield_system(
                K,
                F_x,
                support
            )

            Xfree = LowLevelFEM.solve_linear_system(
                prepared.A,
                prepared.B;
                solver=:backslash
            )

            residual = norm(
                prepared.A * Xfree - prepared.B
            ) / max(norm(prepared.B), eps(Float64))

            @test residual < 1e-9
            @test size(prepared.T, 2) < size(K.A, 1)
            @test all(isfinite, v.a)
            @test all(isfinite, p.a)

            # Stokes/Navier-Stokes block systems are symmetric but indefinite.
            # :auto must therefore fall back from Cholesky to LU.
            vp_sym = @test_logs (
                :warn,
                "Cholesky factorization failed because the symmetric system is not positive definite. Falling back to LU factorization."
            ) solveField(
                Symmetric(K),
                F_x;
                support=support
            )

            v_sym, p_sym = vp_sym

            @test relerr(v_sym.a, v.a) < 1e-9
            @test relerr(p_sym.a, p.a) < 1e-8
        end
    end
end

@testset "Analytical coupled scalar fields" begin
    with_test_gmsh() do
        line_mesh(n=2, order=1)
        mat = Material("body")
        P = Field([mat], type=:ScalarField, dim=1, fieldName=:p)
        Q = Field([mat], type=:ScalarField, dim=1, fieldName=:q)
        A = ∫(P ⋅ 2.0 ⋅ P)
        B = ∫(P ⋅ Q)
        C = ∫(Q ⋅ 3.0 ⋅ Q)
        K = SystemMatrix([A B; B' C])
        F = SystemVector([∫(P ⋅ 4.0), ∫(Q ⋅ 7.0)])
        p, q = solveField(Symmetric(K), F)
        # [2 1; 1 3] * [1; 2] = [4; 7], pointwise on the domain.
        test_constant_field(p, 1.0)
        test_constant_field(q, 2.0)
        @test K.problems == [P, Q]
        @test K.offsets == [0, ndofs(P)]
        @test size(B.A) == (ndofs(P), ndofs(Q))
        @test relerr(K.A * vcat(p.a, q.a), F.a) < 1e-12
        @test Matrix(K.A[1:ndofs(P), ndofs(P)+1:end]) ≈ Matrix(B.A)

        # Reordering whole fields must also reorder the returned fields.
        Kr = SystemMatrix([C B'; B A])
        Fr = SystemVector([∫(Q ⋅ 7.0), ∫(P ⋅ 4.0)])
        qr, pr = solveField(Kr, Fr)
        @test qr.model === Q
        @test pr.model === P
        @test pr.a ≈ p.a atol=1e-12
        @test qr.a ≈ q.a atol=1e-12
    end
end
