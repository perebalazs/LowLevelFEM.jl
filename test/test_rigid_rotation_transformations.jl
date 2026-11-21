# Standalone: julia --project=. -e 'using Test, LowLevelFEM, LinearAlgebra, SparseArrays; include("test/helpers.jl"); include("test/test_rigid_rotation_transformations.jl")'
# Keep this fixture independent of test_transient_remote.jl and its test order.
@testset "Rigid rotation transformations (2D/3D)" begin
    for dim in (2, 3)
        @testset "$(dim)D" begin
            with_test_gmsh() do
                if dim == 2
                    structured_rect_mesh(n=2, rightRemote=true)
                else
                    structured_box_mesh(n=2, rightRemote=true)
                end
                mat = Material("body", E=100.0, ν=0.25)
                materials = [mat, Material("rightremote")]
                U = Field(materials, type=:VectorField, dim=dim, fieldName=:u)
                Φ = Field(materials, type=dim == 2 ? :ScalarField : :VectorField,
                          dim=dim, fieldName=:φ)
                mpc_u = MPC(master="rightremote", slave="right", field=U)
                mpc_φ = MPC(master="rightremote", slave="right", field=Φ)
                R = rigidRotationMap(mpc_u, mpc_φ)
                nu, nφ = ndofs(U), ndofs(Φ)

                # Independent geometry oracle: u_rot = φ × (x - x_master).
                master_tag = LowLevelFEM.getTagForPhysicalName("rightremote")
                masters, xm = gmsh.model.mesh.getNodesForPhysicalGroup(-1, master_tag)
                master = Int(only(masters))
                slave_tag = LowLevelFEM.getTagForPhysicalName("right")
                slaves, xyz = gmsh.model.mesh.getNodesForPhysicalGroup(-1, slave_tag)
                expected_R = zeros(nu, nφ)
                levers = reshape(xyz, 3, :) .- xm
                for (j, node) in enumerate(slaves)
                    rows = dim * (Int(node) - 1) .+ (1:dim)
                    r = levers[:, j]
                    if dim == 2
                        expected_R[rows, Int(node)] = [-r[2], r[1]]
                    else
                        cols = 3 * (Int(node) - 1) .+ (1:3)
                        for k in 1:3
                            axis = zeros(3)
                            axis[k] = 1.0
                            expected_R[rows, cols[k]] = cross(axis, r)
                        end
                    end
                end
                @test R.A ≈ expected_R atol=1e-14
                @test R.model === Φ
                @test R.test_model === U
                @test norm(expected_R) > 0

                # A tangential force on one slave guarantees a nonzero moment.
                j = argmax(vec(sum(abs2, levers; dims=1)))
                node = Int(slaves[j])
                r = levers[:, j]
                force = dim == 2 ? [-r[2], r[1]] : cross([1.0, 2.0, 3.0], r)
                torque = dim == 2 ? [r[1]*force[2] - r[2]*force[1]] : cross(r, force)
                @test norm(torque) > 1e-8
                fu = zeros(nu, 2)
                urows = dim * (node - 1) .+ (1:dim)
                fu[urows, 1] = force
                fu[:, 2] = -2fu[:, 1]
                fφ = zeros(nφ, 2)
                φrows = Φ.pdim * (node - 1) .+ (1:Φ.pdim)
                expected_moment = zeros(nφ, 2)
                expected_moment[φrows, 1] = torque
                expected_moment[:, 2] = -2expected_moment[:, 1]
                @test R.A' * fu ≈ expected_moment atol=1e-13

                Ku = ∫(SymGrad(U) ⋅ D(dim == 2 ? :PlaneStress : :Solid, mat) ⋅ SymGrad(U); Ω="body")
                Z = SystemMatrix(spzeros(nu, nφ), Φ, U)
                Zφ = SystemMatrix(spzeros(nφ, nφ), Φ)
                keys_u = dim == 2 ? (; ux=0.0, uy=0.0) : (; ux=0.0, uy=0.0, uz=0.0)
                support = [BoundaryCondition("left"; field=U, keys_u...)]
                physical_solutions = Matrix{Float64}[]

                for reverse_fields in (false, true)
                    @testset "Field order $(reverse_fields ? "Φ,U" : "U,Φ")" begin
                        K0 = reverse_fields ? SystemMatrix([Zφ Z'; Z Ku]) : SystemMatrix([Ku Z; Z' Zφ])
                        F0 = reverse_fields ? vcat(fφ, fu) : vcat(fu, fφ)
                        f0 = SystemVector(F0, nothing, K0.problems, K0.offsets)
                        original_K, original_F = copy(K0.A), copy(f0.a)
                        iu, iφ = reverse_fields ? (2, 1) : (1, 2)
                        ur = K0.offsets[iu] .+ (1:nu)
                        φr = K0.offsets[iφ] .+ (1:nφ)
                        n = nu + nφ
                        expected_T = Matrix{Float64}(I, n, n)
                        expected_T[ur, φr] = expected_R
                        T = rigidRotationTransformation(K0, R)
                        @test T.T isa SparseMatrixCSC
                        @test T.T ≈ expected_T atol=1e-14
                        @test T.non == n
                        @test T.dim == 1
                        for TT in (transpose(T), adjoint(T))
                            @test TT.T isa SparseMatrixCSC
                            @test TT.T ≈ expected_T'
                            @test (TT.non, TT.dim) == (T.non, T.dim)
                        end

                        K1, f1 = T' * K0 * T, T' * f0
                        expected_f = copy(F0)
                        expected_f[φr, :] += expected_moment
                        @test f1.a ≈ expected_f atol=1e-13
                        @test (T * f0).a ≈ expected_T * F0 atol=1e-13
                        @test size(f1.a) == (n, 2)
                        @test K1.A ≈ expected_T' * K0.A * expected_T
                        @test K1.problems === K0.problems
                        @test K1.offsets === K0.offsets
                        @test f1.model === f0.model
                        @test f1.problems === f0.problems
                        @test f1.offsets === f0.offsets
                        @test K0.A == original_K
                        @test f0.a == original_F

                        q = sin.(collect(1.0:n))
                        physical_q = expected_T * q
                        @test dot(q, K1.A*q) ≈ dot(physical_q, K0.A*physical_q) rtol=1e-12
                        @test vec(q' * f1.a) ≈ vec(physical_q' * F0) atol=1e-12

                        # A force acting only at the master has no R' * fu term.
                        fm = zeros(n, 1)
                        master_rows = dim * (master - 1) .+ (1:dim)
                        fm[ur[master_rows], 1] = force
                        remote = SystemVector(fm, nothing, K0.problems, K0.offsets)
                        @test (T' * remote).a ≈ fm atol=1e-14

                        # Both the matrix shape and the RHS row count are checked.
                        bad_size = Transformation(sparse(I, n-1, n-1), n-1, 1)
                        rectangular = Transformation(sparse(ones(n, n-1)), n, 1)
                        short_rhs = SystemVector(zeros(n-1, 2), nothing, K0.problems, K0.offsets)
                        @test_throws DimensionMismatch bad_size * f0
                        @test_throws DimensionMismatch rectangular * f0
                        @test_throws DimensionMismatch T * short_rhs

                        # solveField applies the MPC after the explicit rigid map.
                        prepared = LowLevelFEM.prepare_multifield_system(
                            K1, f1, support; mpc=[mpc_u, mpc_φ])
                        @test size(prepared.T, 2) < n
                        G = expected_T * prepared.T[:, prepared.free_r]
                        reference_K = G' * K0.A * G
                        reference_f = G' * F0
                        @test prepared.A ≈ reference_K rtol=1e-12 atol=1e-12
                        @test prepared.B ≈ reference_f atol=1e-13
                        @test norm(reference_f - (prepared.T' * F0)[prepared.free_r, :]) > 1e-8
                        solution = solveField(K1, f1; support, mpc=[mpc_u, mpc_φ], solver=:backslash)
                        u, φ = reverse_fields ? reverse(solution) : solution
                        physical = u.a + expected_R * φ.a
                        reference = G * (reference_K \ reference_f)
                        @test physical ≈ reference[ur, :] rtol=1e-9 atol=1e-12
                        for (field, P, mpc) in ((u, U, mpc_u), (φ, Φ, mpc_φ))
                            rep = LowLevelFEM.mpcRepresentativeMap(P, [mpc])
                            tied = findall(i -> rep[i] != i, eachindex(rep))
                            @test !isempty(tied)
                            @test field.a[tied, :] ≈ field.a[rep[tied], :] atol=1e-12
                        end
                        push!(physical_solutions, physical)
                    end
                end
                @test physical_solutions[1] ≈ physical_solutions[2] rtol=1e-9 atol=1e-12
            end
        end
    end
end
