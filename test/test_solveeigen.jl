@testset "solveEigen: analytical FEM eigenproblem" begin
    with_test_gmsh() do
        line_mesh(n=2, order=1)
        P = Field([Material("body")], type=:ScalarField, dim=1, fieldName=:p)
        K = ∫(Grad(P) ⋅ Grad(P))
        M = ∫(P ⋅ P)
        support = [BoundaryCondition("left"; field=P, p=7.0),
                   BoundaryCondition("right"; field=P, p=-3.0)]
        # Two P1 elements: K_free=4, M_free=1/3, hence λ=12.
        modes = solveEigen(Symmetric(K), Symmetric(M); n=1, solver=:direct,
                           support=support, normalize=:mass)
        @test modes isa ScalarField
        @test modes.kind === :eigen_real
        @test eigenValue(modes) ≈ 12.0 atol=1e-10
        @test all(iszero, modes.a[constrainedDoFs(P, support), :])
        @test dot(modes.a[:, 1], M.A * modes.a[:, 1]) ≈ 1.0 atol=1e-12
        free = freeDoFs(P, support)
        test_eigen_residual(K.A, M.A, modes.t, modes.a;
                            transform=sparse(I, ndofs(P), ndofs(P))[:, free])

        selected = eigenVector(modes; n=1)
        @test selected.a == modes.a
        @test selected.kind === modes.kind
        selected.a[free[1], 1] = 99.0
        @test modes.a[free[1], 1] != 99.0

        λ, ϕ = solveEigen(K, M; nev=1, solver=:direct, support=support, output=:raw)
        @test λ ≈ modes.t
        @test norm(ϕ[:, 1]) ≈ 1.0 atol=1e-12
    end
end

@testset "Complex reduced spectrum regression" begin
    # Preserve the complex-valued algebraic oracle from eigen_api_tests.jl.
    # SystemMatrix stores Float64, so this backend-specific case is checked
    # directly at the reduced solver boundary.
    A = spdiagm(0 => ComplexF64[3+2im, -2+1im, 5-3im, 1-1im])
    B = spdiagm(0 => ones(ComplexF64, 4))
    pencil = (; A, B, symmetricA=false, symmetricB=false, uploA=:U, uploB=:U)
    λ, ϕ = LowLevelFEM.solve_eigen_reduced(pencil; n=2, shift=0.0, solver=:direct)
    @test real.(λ) ≈ [1.0, 3.0]
    test_eigen_residual(A, B, λ, ϕ)
    @test_throws ArgumentError LowLevelFEM.solve_eigen_reduced(pencil; n=2,
                                                              shift=0.0, solver=:arpack)
end

@testset "solveEigen: ordered selection and backends" begin
    with_test_gmsh() do
        line_mesh(n=39, order=1)
        P = Field([Material("body")], type=:ScalarField, dim=1, fieldName=:p)
        N = ndofs(P)
        A = SystemMatrix(spdiagm(0 => collect(1.0:N)), P)
        B = SystemMatrix(spdiagm(0 => ones(N)), P)
        # Preserve the exact spectra from examples/dev/eigen/eigen_api_tests.jl.
        for (options, expected) in (((n=4,), [1.0, 2.0, 3.0, 4.0]),
                                    ((n=4, shift=12.0), [13.0, 14.0, 15.0, 16.0]),
                                    ((n=3, sigma=12.2, which=:LM), [11.0, 12.0, 13.0]),
                                    ((n=3, which=:LM), [38.0, 39.0, 40.0]))
            λ, ϕ = solveEigen(A, B; solver=:direct, output=:raw, options...)
            @test λ ≈ expected atol=1e-10
            test_eigen_residual(A.A, B.A, λ, ϕ)
        end
        for (options, expected) in (((n=4, shift=12.0), [13.0, 14.0, 15.0, 16.0]),
                                    ((n=3, sigma=12.2, which=:LM), [11.0, 12.0, 13.0]))
            λ, ϕ = solveEigen(A, B; solver=:arpack, output=:raw, options...,
                              solveroptions=(ncv=24, tol=1e-11, v0=ones(N)))
            @test λ ≈ expected atol=1e-6
            test_eigen_residual(A.A, B.A, λ, ϕ)
        end

        modes = solveEigen(A, B; n=4, solver=:direct)
        @test modes.nsteps == 4
        @test eigenValue(modes; n=3) ≈ 3.0
        @test eigenVector(modes; n=3).a == modes.a[:, 3:3]
        real_modes, imag_modes = solveEigen(A, B; n=4, solver=:direct, part=:both)
        @test real_modes.kind === :eigen_real
        @test imag_modes.kind === :eigen_imag
        @test all(iszero, imag_modes.a)
        @test all(iszero, imag_modes.t)
        @test solveEigen(A, B; n=4, solver=:direct, part=:imag).kind === :eigen_imag

        # Nonidentity B distinguishes B-normalization from Euclidean scaling.
        metric = SystemMatrix(2.0 * B.A, P)
        normalized = solveEigen(A, metric; n=3, solver=:direct, normalize=:B)
        @test normalized.a' * metric.A * normalized.a ≈ Matrix{Float64}(I, 3, 3) atol=1e-12

        # Gmsh views are headless and cleared by with_test_gmsh.
        @test showEigenResults(modes; visible=false) isa Integer
        @test showModalResults(modes; visible=false) isa Integer
        @test showBucklingResults(modes; visible=false) isa Integer
    end
end

@testset "solveEigen: complex pairs and multifield output" begin
    with_test_gmsh() do
        structured_rect_mesh(n=1, order=1)
        mat = Material("body")
        P = Field([mat], type=:VectorField, dim=2, fieldName=:u)
        N = ndofs(P)
        A = spdiagm(0 => collect(2.0:N+1.0))
        A[1, 1] = A[2, 2] = 0.0
        A[1, 2], A[2, 1] = -1.0, 1.0
        K, M = SystemMatrix(A, P), SystemMatrix(spdiagm(0 => ones(N)), P)
        λ, ϕ = solveEigen(K, M; n=2, solver=:direct, output=:raw)
        @test sort(imag.(λ)) ≈ [-1.0, 1.0] atol=1e-12
        test_eigen_residual(A, M.A, λ, ϕ)
        real_modes, imag_modes = solveEigen(K, M; n=2, solver=:direct, part=:both)
        @test real_modes isa VectorField
        @test imag_modes isa VectorField
        @test real_modes.type === :v2D
        @test real_modes.kind === :eigen_real
        @test imag_modes.kind === :eigen_imag
        test_eigen_residual(A, M.A, real_modes.t + im * imag_modes.t,
                            real_modes.a + im * imag_modes.a)

        S = Field([mat], type=:ScalarField, dim=2, fieldName=:s)
        Ks = ∫(S ⋅ 2.0 ⋅ S)
        Kv = ∫(P ⋅ 5.0 ⋅ P)
        Ms, Mv = ∫(S ⋅ S), ∫(P ⋅ P)
        Zsv = SystemMatrix(spzeros(ndofs(S), N), P, S)
        Zvs = Zsv'
        Kblock = SystemMatrix([Ks Zsv; Zvs Kv])
        Mblock = SystemMatrix([Ms Zsv; Zvs Mv])
        s, v = solveEigen(Kblock, Mblock; n=ndofs(S)+1, solver=:direct)
        @test s isa ScalarField
        @test v isa VectorField
        @test s.t == v.t
        @test s.t ≈ vcat(fill(2.0, ndofs(S)), 5.0) atol=1e-10
        @test size(s.a) == (ndofs(S), ndofs(S)+1)
        @test size(v.a) == (N, ndofs(S)+1)
        test_eigen_residual(Kblock.A, Mblock.A, s.t, vcat(s.a, v.a))
        sr, si, vr, vi = solveEigen(Kblock, Mblock; n=2, solver=:direct, part=:both)
        @test (sr.model, si.model, vr.model, vi.model) == (S, S, P, P)
        @test all(iszero, si.a) && all(iszero, vi.a)

        support = [BoundaryCondition("left"; field=S, s=1.0),
                   BoundaryCondition("right"; field=P, ux=0.0, uy=0.0)]
        s_bc, v_bc = solveEigen(Kblock, Mblock; n=2, solver=:direct, support=support)
        @test all(iszero, s_bc.a[constrainedDoFs(S, support), :])
        @test all(iszero, v_bc.a[constrainedDoFs(P, support), :])
    end
end

@testset "solveEigen: MPC, coordinates and reduced interpolation" begin
    with_test_gmsh() do
        line_mesh(n=4, order=1)
        P = Field([Material("body")], type=:ScalarField, dim=1, fieldName=:p)
        K, M = ∫(Grad(P) ⋅ Grad(P)), ∫(P ⋅ P)
        mpc = MPC(master="left", slave="right", field=P, p=true)
        modes = solveEigen(K, M; n=3, solver=:direct, mpc=[mpc])
        left = constrainedDoFs(P, [BoundaryCondition("left"; field=P, p=0.0)])
        right = constrainedDoFs(P, [BoundaryCondition("right"; field=P, p=0.0)])
        @test modes.a[left, :] ≈ modes.a[right, :] atol=1e-12
        @test abs(first(modes.t)) < 1e-10
        prepared = LowLevelFEM.prepare_eigen_system(K, M; mpc=[mpc])
        test_eigen_residual(K.A, M.A, modes.t[2:end], modes.a[:, 2:end];
                            transform=prepared.T[:, prepared.free])
    end

    with_test_gmsh() do
        structured_rect_mesh(n=1, order=2)
        mat = Material("body")
        U = Field([mat], type=:VectorField, dim=2, fieldName=:u)
        K, M = ∫(SymGrad(U) ⋅ D(:PlaneStress, mat) ⋅ SymGrad(U)), ∫(U ⋅ U)
        support = [BoundaryCondition("left"; field=U, ux=0.0, uy=0.0)]
        e1 = VectorField(U, "body", [0.0, 1.0, 0.0])
        cs = CoordinateSystem(e1)
        reference = solveEigen(K, M; n=3, solver=:direct, support=support)
        rotated = solveEigen(K, M; n=3, solver=:direct, support=support, coordSys=[cs])
        @test rotated.t ≈ reference.t rtol=1e-10

        R = Field([mat], type=:VectorField, dim=2, fieldName=:r, reducedOrder=true)
        Kr, Mr = ∫(SymGrad(R) ⋅ D(:PlaneStress, mat) ⋅ SymGrad(R)), ∫(R ⋅ R)
        bc = [BoundaryCondition("left"; field=R, rx=0.0, ry=0.0)]
        modes = solveEigen(Kr, Mr; n=2, solver=:direct, support=bc)
        T, restriction = reductionMatrices(R)
        @test size(T, 2) < size(T, 1)
        @test modes.a ≈ T * (restriction * modes.a) atol=1e-10
        @test all(iszero, modes.a[constrainedDoFs(R, bc), :])
        prepared = LowLevelFEM.prepare_eigen_system(Kr, Mr; support=bc)
        test_eigen_residual(Kr.A, Mr.A, modes.t, modes.a;
                            transform=prepared.T[:, prepared.free])
    end
end
