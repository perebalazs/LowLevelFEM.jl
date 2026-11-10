@testset "Preconditioners respect the authoritative symmetric triangle" begin
    A = sparse([4.0 1.0; 1.0 3.0])
    b = [1.0, 2.0]
    for wrapper in (Symmetric, Hermitian), uplo in (:U, :L)
        storage = copy(A)
        storage[uplo === :U ? CartesianIndex(2,1) : CartesianIndex(1,2)] = 99.0
        S = wrapper(storage, uplo)
        for pc in (:ilu, :ich)
            P = LowLevelFEM._build_preconditioner(Val(pc), S)
            @test P \ b ≈ A \ b atol=1e-12
            @test LowLevelFEM.solve_linear_system(S, b; solver=:cg,
                preconditioner=P, solveroptions=(reltol=1e-12,)) ≈ A \ b
        end
    end
    # Hermitian materialization must conjugate the ignored triangle.
    H = Hermitian(sparse([4.0+0im 1.0+1im; 99.0+0im 3.0+0im]), :U)
    P = LowLevelFEM._build_preconditioner(Val(:ilu), H)
    @test P \ complex.(b) ≈ H \ b atol=1e-12
end

@testset "Iterative multiple RHS use a solution-compatible element type" begin
    A = sparse([4.0 1.0; 1.0 3.0])
    B = [1 2; 2 1]
    for solver in (:cg, :gmres), rhs in (B, sparse(B))
        prepared = LowLevelFEM.prepare_linear_solver(A; solver,
            preconditioner=:none, solveroptions=(reltol=1e-12,))
        @test LowLevelFEM.solve_prepared_linear_system(prepared, rhs) ≈ A \ B
        @test size(LowLevelFEM.solve_prepared_linear_system(prepared, zeros(Int,2,0))) == (2,0)
    end
    for T in (Float32, Float64, ComplexF64)
        a, b = T.(A), T.(B)
        tolerance = typeof(real(zero(T)))(1e-6)
        x = LowLevelFEM.solve_linear_system(a, b; solver=:cg,
            preconditioner=:none, solveroptions=(reltol=tolerance,))
        @test x ≈ a \ b rtol=1e-5
    end
end

@testset "Iterative nonconvergence is explicit" begin
    A, b = sparse([4.0 1.0; 1.0 3.0]), [1.0, 2.0]
    for solver in (:cg, :gmres)
        options = (; maxiter=1, reltol=1e-14)
        failure = try
            LowLevelFEM.solve_linear_system(A, b; solver,
                preconditioner=:none, solveroptions=options)
            nothing
        catch err
            err
        end
        @test failure isa LinearSolveConvergenceError
        if failure isa LinearSolveConvergenceError
            @test failure.solver === solver
            @test failure.iterations == 1
            @test failure.residual_norm > 1e-3
            @test occursin("did not converge", sprint(showerror, failure))
        end
        @test_logs (:warn, r"did not converge") begin
            x = LowLevelFEM.solve_linear_system(A, b; solver, preconditioner=:none,
                solveroptions=merge(options, (; on_failure=:warn)))
            @test norm(A*x-b) > 1e-3
        end
        @test_throws ArgumentError LowLevelFEM.prepare_linear_solver(A; solver,
            preconditioner=:none, solveroptions=(on_failure=:ignore,))
    end

    with_test_gmsh() do
        line_mesh(n=1, order=1)
        P = Field([Material("body")], type=:ScalarField, dim=1, fieldName=:p)
        K = SystemMatrix(sparse([4.0 1.0; 1.0 3.0]), P)
        M = SystemMatrix(spdiagm(0=>ones(2)), P)
        load = ScalarField(Matrix{Float64}[], reshape(b,:,1), [0.0], Int[], 1, :scalar, P)
        zero = ScalarField(Matrix{Float64}[], zeros(2,1), [0.0], Int[], 1, :scalar, P)
        options = (; maxiter=1, reltol=1e-14)
        @test_throws LinearSolveConvergenceError solveField(K, load; solver=:cg,
            preconditioner=:none, solveroptions=options)
        @test_throws LinearSolveConvergenceError FDM(K, M, load, zero, 3, 0.1;
            solver=:cg, preconditioner=:none, solveroptions=options)
    end
end
