@testset "Single-field boundary conditions and MPC" begin
    with_test_gmsh() do
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
    end
end

@testset "Multifield boundary conditions and MPC" begin
    with_test_gmsh() do
        (; P, V, A, B, C, D, K, fv_x, fp, F_x, pres, supp_top, supp_bottom, support) = stokes_fixture()
        mpc_v = MPC(
            master="rightbottom",
            slave="right",
            field=V,
            vx=true,
            vy=true
        )

        mpc_p = MPC(
            master="rightbottom",
            slave="leftbottom",
            field=P,
            p=true
        )

        mpcs = [mpc_v, mpc_p]

        @testset "Multifield MPC with reduced-order pressure" begin
            v_mpc, p_mpc = solveField(
                K,
                F_x;
                support=support,
                mpc=mpcs
            )

            rep_v = LowLevelFEM.mpcRepresentativeMap(V, [mpc_v])
            rep_p = LowLevelFEM.mpcRepresentativeMap(P, [mpc_p])

            tied_v = [i for i in eachindex(rep_v) if rep_v[i] != i]
            tied_p = [i for i in eachindex(rep_p) if rep_p[i] != i]

            @test !isempty(tied_v)
            @test !isempty(tied_p)

            @test maximum(
                abs(v_mpc.a[i, 1] - v_mpc.a[rep_v[i], 1])
                for i in tied_v
            ) < 1e-10

            @test maximum(
                abs(p_mpc.a[i, 1] - p_mpc.a[rep_p[i], 1])
                for i in tied_p
            ) < 1e-10
        end

        # 90-degree local basis on the complete velocity field:
        # local e1 = global +y, local e2 = global -x.
        e1_v = VectorField(
            V,
            "body",
            [0.0, 1.0, 0.0]
        )

        cs_v = CoordinateSystem(e1_v)

        fv_global_y = ∫(V ⋅ [0.0, 1.0])

        F_global_y = SystemVector([fv_global_y, fp])

        fv_local_e1 = ∫(V ⋅ [1.0, 0.0])

        F_local_e1 = SystemVector([fv_local_e1, fp])

        @testset "Multifield CoordinateSystem and Global(SystemVector)" begin
            v_ref, p_ref = solveField(
                K,
                F_global_y;
                support=support
            )

            v_local, p_local = solveField(
                K,
                F_local_e1;
                support=support,
                coordSys=[cs_v]
            )

            v_global, p_global = solveField(
                K,
                Global(F_global_y);
                support=support,
                coordSys=[cs_v]
            )

            @test relerr(v_local.a, v_ref.a) < 1e-9
            @test relerr(p_local.a, p_ref.a) < 1e-8
            @test relerr(v_global.a, v_ref.a) < 1e-9
            @test relerr(p_global.a, p_ref.a) < 1e-8

            prepared = LowLevelFEM.prepare_multifield_system(
                K,
                F_local_e1,
                support;
                coordSys=[cs_v]
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
            @test P.reducedOrder
            @test prepared.Q !== nothing
        end

        @testset "CoordinateSystem + MPC" begin
            v_ref, p_ref = solveField(
                K,
                F_global_y;
                support=support,
                mpc=mpcs
            )

            v_qmpc, p_qmpc = solveField(
                K,
                F_local_e1;
                support=support,
                mpc=mpcs,
                coordSys=[cs_v]
            )

            @test relerr(v_qmpc.a, v_ref.a) < 1e-9
            @test relerr(p_qmpc.a, p_ref.a) < 1e-8

            # Also exercise the single-field path with the same velocity field.
            Kv = A + D

            v_single_ref = solveField(
                Kv,
                fv_global_y;
                support=[supp_top, supp_bottom],
                mpc=[mpc_v]
            )

            v_single_qmpc = solveField(
                Kv,
                fv_local_e1;
                support=[supp_top, supp_bottom],
                mpc=[mpc_v],
                coordSys=[cs_v]
            )

            @test relerr(
                v_single_qmpc.a,
                v_single_ref.a
            ) < 1e-9
        end

        @testset "Rotated MPC / sector-symmetry relation" begin
            θ = π / 6
            Lx = 2.0

            sector_e1(x, y, z) = begin
                α = θ * x / Lx
                [cos(α), sin(α), 0.0]
            end

            e1_sector = VectorField(
                V,
                "body",
                sector_e1
            )

            cs_sector = CoordinateSystem(e1_sector)

            mpc_sector = MPC(
                master="leftbottom",
                slave="right",
                field=V,
                vx=true,
                vy=true
            )

            support_sector = [
                BoundaryCondition("left", field=V, vx=0.0),
                BoundaryCondition("lefttop", field=V, vy=0.0),
                pres
            ]

            v_sector, _ = solveField(
                K,
                Global(F_global_y);
                support=support_sector,
                mpc=[mpc_sector],
                coordSys=[cs_sector]
            )

            prepared = LowLevelFEM.prepare_multifield_system(
                K,
                Global(F_global_y),
                support_sector;
                mpc=[mpc_sector],
                coordSys=[cs_sector]
            )

            nV = V.non * V.pdim
            Qv = prepared.Q[1:nV, 1:nV]
            v_local = Qv' * v_sector.a[:, 1]

            rep = LowLevelFEM.mpcRepresentativeMap(
                V,
                [mpc_sector]
            )

            tied = [i for i in eachindex(rep) if rep[i] != i]
            @test !isempty(tied)

            local_error = maximum(
                abs(v_local[i] - v_local[rep[i]])
                for i in tied
            )

            @test local_error < 1e-10

            rotation_errors = Float64[]

            for node in 1:V.non
                sdofs = [2node - 1, 2node]
                mdofs = [rep[sdofs[1]], rep[sdofs[2]]]

                if mdofs != sdofs
                    us = v_sector.a[sdofs, 1]
                    um = v_sector.a[mdofs, 1]
                    Qs = Matrix(Qv[sdofs, sdofs])
                    Qm = Matrix(Qv[mdofs, mdofs])

                    push!(
                        rotation_errors,
                        norm(us - Qs * Qm' * um)
                    )
                end
            end

            @test !isempty(rotation_errors)
            @test maximum(rotation_errors) < 1e-10
        end
    end
end
