# Regression tests for multifield solveField transformations.

@testset "solveField multifield" begin
    initialized_here = gmsh.isInitialized() == 0
    initialized_here && gmsh.initialize()

    try
        relerr(a, b) = norm(a - b) / max(norm(b), eps(Float64))

        structured_rect_mesh(
            lx=2.0,
            n=2,
            order=2
        )

        material = Material("body")

        P = Field(
            [material],
            type=:ScalarField,
            fieldName=:p,
            rhsName=:fp,
            dim=2,
            reducedOrder=true
        )

        V = Field(
            [material],
            type=:VectorField,
            fieldName=:v,
            rhsName=:fv,
            dim=2
        )

        pres = BoundaryCondition(
            "rightbottom",
            field=P,
            p=0.0
        )

        supp_top = BoundaryCondition(
            "top",
            field=V,
            vx=0.0,
            vy=0.0
        )

        supp_bottom = BoundaryCondition(
            "bottom",
            field=V,
            vx=0.0,
            vy=0.0
        )

        support = [
            pres,
            supp_top,
            supp_bottom
        ]

        fv_x = loadVector(
            V,
            [LoadCondition("body", fvx=1.0, fvy=0.0)]
        )

        fp = loadVector(P, [])
        F_x = SystemVector([fv_x, fp])

        μ = 1.0
        γ = 1e-1

        A = ∫((SymGrad(V) ⋅ SymGrad(V)) * 2μ)
        B = ∫(Div(V) ⋅ P)
        C = 0 * ∫(P ⋅ P)
        D = ∫(Div(V) ⋅ Div(V) * γ)

        K = SystemMatrix([
            A + D  B
            B'     C
        ])

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
                :info,
                "Symmetric system is not positive definite; using LU factorization."
            ) solveField(
                Symmetric(K),
                F_x;
                support=support
            )

            v_sym, p_sym = vp_sym

            @test relerr(v_sym.a, v.a) < 1e-9
            @test relerr(p_sym.a, p.a) < 1e-8
        end

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

        fv_global_y = loadVector(
            V,
            [LoadCondition("body", fvx=0.0, fvy=1.0)]
        )

        F_global_y = SystemVector([fv_global_y, fp])

        fv_local_e1 = loadVector(
            V,
            [LoadCondition("body", fvx=1.0, fvy=0.0)]
        )

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

    finally
        initialized_here && gmsh.finalize()
    end
end
