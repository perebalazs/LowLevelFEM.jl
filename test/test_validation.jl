@testset "Validation and error handling" begin
    with_test_gmsh() do
        structured_rect_mesh(n=1, order=1)
        mat = Material("body")
        P = Field([mat], type=:ScalarField, dim=2, fieldName=:p)
        Q = Field([mat], type=:ScalarField, dim=2, fieldName=:q)
        K = ∫(P ⋅ P)
        f = ∫(P ⋅ 1.0)

        @testset "Assembly inputs" begin
            @test_throws ErrorException ∫(Grad(P) ⋅ ones(3, 2) ⋅ Grad(P))
            @test_throws ErrorException ∫(P ⋅ P; Ω="body", Γ="right")
            @test_throws ErrorException ∫(P ⋅ P; Ω="right")
            @test_throws ErrorException ∫(P ⋅ P; Γ="body")
            @test_throws ErrorException ∫(P ⋅ P; Ω="missing_group")
            @test_throws ErrorException ∫(P ⋅ P; assembly=:unknown)
            for threads in (0, true, :unknown, Threads.nthreads()+1)
                @test_throws ArgumentError ∫(P ⋅ P; threads=threads)
            end
            pattern = build_csc_pattern(P, P)
            @test_throws ErrorException ∫(P ⋅ P; assembly=:ijv, csc_matrix=pattern)
            @test_throws ErrorException ∫(P ⋅ P; csc_matrix=spzeros(1, 1))
        end

        @testset "Field and RHS compatibility" begin
            nodal = scalarField(P, "body", 1.0)
            elementwise = ScalarField(P, "body", 1.0)
            @test_throws ErrorException nodal + elementwise
            @test_throws ArgumentError ScalarField([], zeros(ndofs(P), 1), [0.0],
                                                  Int[], 1, :scalar, P, :invalid)
            @test_throws ErrorException solveField(K, ∫(Q ⋅ 1.0))
            short = ScalarField([], zeros(ndofs(P)-1, 1), [0.0], Int[], 1, :scalar, P)
            @test_throws ErrorException solveField(K, short)
            @test_throws ErrorException solveField(K, f; solver=:unknown)
            @test_throws ErrorException solveField(K, f; solveroptions=(unsupported=1,))
        end

        @testset "Eigen choices, metadata and accessors" begin
            for options in ((n=0,), (n=ndofs(P)+1,), (n=1, nev=1),
                            (solver=:unknown,), (which=:UNKNOWN,),
                            (shift=1.0, sigma=1.0), (shift=1.0, which=:LM),
                            (shift=Inf,), (normalize=:unknown,),
                            (output=:unknown,), (part=:unknown,),
                            (solveroptions=(invalid=4,),),
                            (ritzvec=false,), (unsupported=1,),
                            (tol=1e-8, solveroptions=(tol=1e-9,)),
                            (solveroptions=[1, 2],))
                @test_throws ArgumentError solveEigen(K, K; n=1, options...)
            end
            @test_throws ErrorException solveEigen(K, ∫(Q ⋅ Q); n=1, solver=:direct)
            small = SystemMatrix(spdiagm(0 => ones(ndofs(P)-1)), P)
            @test_throws ErrorException solveEigen(K, small; n=1, solver=:direct)
            support = [BoundaryCondition("body"; field=P, p=0.0)]
            @test_throws ErrorException solveEigen(K, K; n=1, support=support, solver=:direct)
            @test_throws ArgumentError eigenValue(f)
            @test_throws ArgumentError eigenVector(f)
            modes = solveEigen(K, K; n=1, solver=:direct)
            @test_throws ArgumentError eigenValue(modes; n=0)
            @test_throws ArgumentError eigenVector(modes; n=2)
            @test_throws ArgumentError mergeFields([modes, modes])
            @test_throws ArgumentError eigenValue(modes + modes)
        end

        @testset "Multifield metadata ordering" begin
            Mq = ∫(Q ⋅ Q)
            B = ∫(P ⋅ Q)
            block = SystemMatrix([K B; B' Mq])
            wrong_order = SystemVector([∫(Q ⋅ 1.0), f])
            @test_throws ErrorException solveField(block, wrong_order)
            wrong = SystemMatrix(spzeros(1, 1), Q, P)
            @test_throws ErrorException SystemMatrix([K wrong; wrong' Mq])
        end

        @testset "Conflicting MPC prescribed values" begin
            mpc = MPC(master="left", slave="right", field=P, p=true)
            support = [BoundaryCondition("left"; field=P, p=0.0),
                       BoundaryCondition("right"; field=P, p=1.0)]
            @test_throws ErrorException solveField(K, f; support=support, mpc=[mpc])
        end

        @testset "Coordinate basis validation" begin
            U = Field([mat], type=:VectorField, dim=2, fieldName=:u)
            Ku, fu = ∫(U ⋅ U), ∫(U ⋅ [1.0, 0.0])
            zero_direction = CoordinateSystem(VectorField(U, "body", [0.0, 0.0, 0.0]))
            @test_throws ErrorException solveField(Ku, fu; coordSys=[zero_direction])
        end
    end

    with_test_gmsh() do
        structured_box_mesh(n=1, order=2)
        mat = Material("body")
        U = Field([mat], type=:VectorField, dim=3, fieldName=:u, reducedOrder=true)
        K = ∫(SymGrad(U) ⋅ D(:Solid, mat) ⋅ SymGrad(U))
        M = ∫(U ⋅ U)
        f = ∫(U ⋅ [1.0, 0.0, 0.0]; Γ="right")
        support = [BoundaryCondition("left"; field=U, ux=0.01, uy=0.0, uz=0.0)]
        cs = CoordinateSystem(VectorField(U, "right", [0.0, 1.0, 0.0]),
                              VectorField(U, "right", [-1.0, 0.0, 0.0]))
        # Preserve the original explicit unsupported-combination regression.
        err = try
            solveField(K, f; support=support, coordSys=[cs])
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("not yet implemented", sprint(showerror, err))
        @test_throws ErrorException solveEigen(K, M; n=1, support=support,
                                               coordSys=[cs], solver=:direct)
    end

    with_test_gmsh() do
        (; U, C, zero_u) = flat_contact_fixture()
        @test_throws ErrorException ContactGap(C; components=:invalid)
        @test_throws ErrorException ContactGap(C; active=:invalid)
        for options in ((normal_sign=0.0,), (leaf_size=0,), (projection_tol=0.0,),
                        (projection_maxiter=0,), (aabb_padding=-1.0,),
                        (self_exclusion_layers=-1,), (topology_angle=181.0,))
            @test_throws ErrorException contact(U; master="master", slave="slave", options...)
        end
    end
end
