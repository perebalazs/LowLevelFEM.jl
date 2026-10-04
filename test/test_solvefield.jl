# Regression tests for the unified solveField pipeline.

@testset "solveField single-field" begin
    initialized_here = gmsh.isInitialized() == 0
    initialized_here && gmsh.initialize()

    try
        relerr(a, b) = norm(a - b) / max(norm(b), eps(Float64))

        @testset "RHS handling and linear solvers" begin
            using IterativeSolvers, Preconditioners, IncompleteLU
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

            # Backward-compatible keywords.
            u_legacy_iterative = solveField(
                Symmetric(K),
                f1;
                support=support,
                iterative=true,
                solveroptions=(reltol=1e-10,maxiter=size(K.A, 1))
            )

            u_legacy_ordering = solveField(
                K,
                f1;
                support=support,
                ordering=false
            )

            @test relerr(u_legacy_iterative.a, u_cg.a) < 1e-12
            @test relerr(u_legacy_ordering.a, u1.a) < 1e-10
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

            # CoordinateSystem + reducedOrder on the same field is deliberately
            # deferred and must fail with an explicit message.
            e1 = VectorField(U, "right", [0.0, 1.0, 0.0])
            e2 = VectorField(U, "right", [-1.0, 0.0, 0.0])
            cs = CoordinateSystem(e1, e2)

            err = try
                solveField(
                    K,
                    f;
                    support=support,
                    coordSys=[cs]
                )
                nothing
            catch e
                e
            end

            @test err isa ErrorException
            @test occursin(
                "not yet implemented",
                sprint(showerror, err)
            )
        end

        @testset "MPC remote point" begin
            structured_rect_mesh(n=2)

            remote_tag = gmsh.model.occ.addPoint(1.2, 0.5, 0.0)
            gmsh.model.occ.synchronize()
            gmsh.model.addPhysicalGroup(0, [remote_tag], -1, "remote")

            gmsh.model.mesh.clear()
            gmsh.model.mesh.generate(2)

            mat = Material("body")
            remote = Material("remote")

            U = Field(
                [mat, remote],
                type=:VectorField,
                dim=2,
                fieldName=:u,
                rhsName=:f
            )

            C = [
                2.0 1.0 0.0
                1.0 2.0 0.0
                0.0 0.0 1.0
            ]

            K = ∫(SymGrad(U) ⋅ C ⋅ SymGrad(U), Ω="body")
            f = ∫(U ⋅ [1.0, 0.0], Γ="right")

            mpc = MPC(
                master="remote",
                slave="right",
                field=U,
                ux=true,
                uy=false
            )

            support = [
                BoundaryCondition("left", ux=0.0, uy=0.0),
                BoundaryCondition("remote", uy=0.0)
            ]

            u = solveField(
                K,
                f;
                support=support,
                mpc=[mpc]
            )

            rep = LowLevelFEM.mpcRepresentativeMap(U, [mpc])
            tied = [i for i in eachindex(rep) if rep[i] != i]

            @test !isempty(tied)
            @test maximum(
                abs(u.a[i, 1] - u.a[rep[i], 1])
                for i in tied
            ) < 1e-10

            # Non-homogeneous prescribed master value must propagate to slaves.
            support_nonhom = [
                BoundaryCondition("left", ux=0.0, uy=0.0),
                BoundaryCondition("remote", ux=0.01, uy=0.0)
            ]

            u_nh = solveField(
                K,
                0.0 * f;
                support=support_nonhom,
                mpc=[mpc]
            )

            @test maximum(
                abs(u_nh.a[i, 1] - u_nh.a[rep[i], 1])
                for i in tied
            ) < 1e-10

            remote_ux = constrainedDoFs(
                U,
                [BoundaryCondition("remote", ux=0.0)]
            )

            @test length(remote_ux) == 1
            @test isapprox(
                u_nh.a[only(remote_ux), 1],
                0.01;
                atol=1e-12,
                rtol=0
            )
        end

        @testset "Coordinate systems and Global RHS" begin
            structured_box_mesh(n=2)

            mat = Material("body")
            U = Field(
                [mat],
                type=:VectorField,
                dim=3,
                fieldName=:u
            )

            K = ∫(ε(U) ⋅ D(:Solid, mat) ⋅ ε(U))

            support = [
                BoundaryCondition("left", ux=0.0, uy=0.0, uz=0.0)
            ]

            # 90-degree basis on the right face:
            # local e1 = global +y, local e2 = global -x.
            e1 = VectorField(U, "right", [0.0, 1.0, 0.0])
            e2 = VectorField(U, "right", [-1.0, 0.0, 0.0])
            cs = CoordinateSystem(e1, e2)

            f_local = ∫(U ⋅ [1.0, 0.0, 0.0], Γ="right")
            f_global = ∫(U ⋅ [0.0, 1.0, 0.0], Γ="right")

            u_local = solveField(
                K,
                f_local;
                support=support,
                coordSys=[cs]
            )

            u_global = solveField(
                K,
                Global(f_global);
                support=support,
                coordSys=[cs]
            )

            @test relerr(u_local.a, u_global.a) < 1e-12

            # Mixed local/global RHS must obey superposition.
            f_body_global = ∫(
                U ⋅ [0.0, 0.0, -1.0],
                Ω="body"
            )

            u_body = solveField(
                K,
                Global(f_body_global);
                support=support,
                coordSys=[cs]
            )

            u_mixed = solveField(
                K,
                f_local + Global(f_body_global);
                support=support,
                coordSys=[cs]
            )

            @test relerr(
                u_mixed.a,
                u_local.a + u_body.a
            ) < 1e-10

            # Inclined roller: local ux=0 on a 45-degree basis.
            e1_45 = VectorField(U, "right", [1.0, 1.0, 0.0])
            e2_45 = VectorField(U, "right", [-1.0, 1.0, 0.0])
            cs45 = CoordinateSystem(e1_45, e2_45)

            support45 = [
                BoundaryCondition("left", ux=0.0, uy=0.0, uz=0.0),
                BoundaryCondition("right", ux=0.0)
            ]

            u45 = solveField(
                K,
                Global(f_global);
                support=support45,
                coordSys=[cs45]
            )

            right_nodes = LowLevelFEM._nodesOnPhysicalGroup(U, "right")
            ux_dofs = 3 .* (right_nodes .- 1) .+ 1
            uy_dofs = 3 .* (right_nodes .- 1) .+ 2

            roller_error = maximum(
                abs,
                u45.a[ux_dofs, 1] + u45.a[uy_dofs, 1]
            )

            scale = maximum(abs, u45.a[:, 1])

            @test roller_error / max(scale, eps(Float64)) < 1e-10

            # A normal-vector based coordinate system must be orthogonal.
            en = normalVector(U, "right")
            et = VectorField(U, "right", [0.0, 1.0, 0.0])
            cs_normal = CoordinateSystem(en, et)

            Q = LowLevelFEM._build_coordinate_transformation(
                U,
                [cs_normal]
            )

            Iq = spdiagm(
                0 => ones(Float64, size(Q.T, 1))
            )

            @test norm(Q.T' * Q.T - Iq) / norm(Iq) < 1e-12
        end

    finally
        initialized_here && gmsh.finalize()
    end
end
