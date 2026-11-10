@testset "Coefficient-chain numerical compatibility" begin
    L = LowLevelFEM
    A = [1.2 0.3; -0.1 1.5]
    B = [0.9 -0.2; 0.4 1.1]
    scalar = Dict(1 => [1.0; 2.0;;], 2 => [3.0; 4.0;;])
    variable = Any[scalar 0.2; -0.1 scalar]
    rectangular = [1.0 0.2 0.3; 0.4 1.1 0.5]
    chains = (
        Any[A, B, A], Any[2.0, A, B, scalar],
        Any[scalar, A, B], Any[A, scalar, B], Any[A, B, scalar],
        Any[variable, A, scalar, B, variable],
        Any[scalar, rectangular, rectangular', scalar],
        Any[scalar, 2.0, scalar],
        Any[isodd(i) ? A/2 : scalar for i in 1:25],
    )
    ws = L.OperatorWorkspace(2, 2, 1, 2, 2, 2, 2)
    for chain in chains
        saved = deepcopy(chain)
        plan = L._prepare_coefficient_chain(chain)
        for elem in (1, 2), h in ([0.25, 0.75], [0.6, 0.4])
            reference = L._eval_coefficient_chain_at_gp(chain, elem, h)
            result = plan isa AbstractVector ?
                L._eval_coefficient_chain_at_gp!(ws, plan, elem, h) :
                L._eval_coefficient_at_gp(plan, elem, h)
            @test result ≈ reference rtol=2e-12 atol=1e-12
        end
        @test chain == saved
    end

    @testset "Preserve range across variable factors" begin
        for (scale, factor) in ((1e-200, 1e200), (1e200, 1e-200))
            chain = Any[Dict(1 => fill(scale, 1, 1)),
                        Matrix(Diagonal(fill(factor, 2))),
                        Matrix(Diagonal(fill(factor, 2)))]
            reference = L._eval_coefficient_chain_at_gp(chain, 1, [1.0])
            plan = L._prepare_coefficient_chain(chain)
            result = L._eval_coefficient_chain_at_gp!(ws, plan, 1, [1.0])
            @test all(isfinite, result)
            @test all(!iszero, diag(result))
            @test result == reference
        end
        for value in (Inf, NaN)
            chain = Any[[value 0.0; 0.0 1.0], A]
            plan = L._prepare_coefficient_chain(chain)
            reference = L._eval_coefficient_chain_at_gp(chain, 1, [1.0])
            @test isequal(L._eval_coefficient_chain_at_gp!(ws, plan, 1, [1.0]), reference)
        end
    end
end

@testset "Bilinear chains preserve fields, loads and quadrature" begin
    with_test_gmsh() do
        structured_rect_mesh(n=2, order=1)
        P = Field([Material("body")], type=:ScalarField, dim=2)
        U = Field([Material("body")], type=:VectorField, dim=2)
        c = ScalarField(P, "body", (x,y,z)->1+x+y)
        s = scalarField(P, "body", (x,y,z)->2+x+2y)
        A = [2.0 0.3; -0.2 1.5]
        B = [1.0 -0.1; 0.4 1.2]
        D = A * B
        reference_matrix = [D[1,1]*c D[1,2]*c; D[2,1]*c D[2,2]*c]
        forms = (Grad(P) ⋅ c ⋅ A ⋅ B ⋅ Grad(P),
                 Grad(P) ⋅ A ⋅ c ⋅ B ⋅ Grad(P),
                 Grad(P) ⋅ A ⋅ B ⋅ c ⋅ Grad(P))
        saved = (copy(A), copy(B), deepcopy(c.A), copy(s.a))
        for gauss in (:full, :reduced, 2), workers in unique((1, Threads.nthreads()))
            reference = ∫(Grad(P) ⋅ reference_matrix ⋅ Grad(P);
                          gauss, threads=1, weight=s).A
            for form in forms
                @test ∫(form; gauss, threads=workers, weight=s).A ≈ reference rtol=2e-12 atol=1e-12
                pattern = build_csc_pattern(P, P)
                fill!(pattern.nzval, 0.125)
                seed = copy(pattern)
                @test ∫(form; gauss, threads=workers, weight=s, csc_matrix=pattern).A ≈
                    seed + reference rtol=2e-12 atol=1e-12
            end
        end
        @test (A, B, c.A, s.a) == saved

        form = first(forms)
        before = copy(∫(form).A)
        c.A .*= 2
        @test ∫(form).A ≈ 2before rtol=2e-12 atol=1e-12

        # The shared coefficient preparer also handles RHS component vectors.
        # They must retain their components instead of becoming scalar products.
        load = ∫(U ⋅ [1.0, 2.0]).a
        @test sum(load[1:2:end]) ≈ 1.0
        @test sum(load[2:2:end]) ≈ 2.0
        variable_load = ∫(U ⋅ [s, 2.0]).a
        @test sum(variable_load[1:2:end]) ≈ 3.5
        @test sum(variable_load[2:2:end]) ≈ 2.0
    end
end
