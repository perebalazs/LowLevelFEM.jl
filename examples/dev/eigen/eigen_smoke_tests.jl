# Run from the package environment after copying the refactored source files:
#   julia --project=. eigen_smoke_tests.jl
# These tests need a working Gmsh installation. They use small algebraic FE
# matrices intentionally, so assembly methods are not part of the test.
using Test, LinearAlgebra, SparseArrays
using LowLevelFEM

function check_residual(A, B, lambda, phi, free; tol=1e-8)
    for j in eachindex(lambda)
        r = A * phi[:, j] - lambda[j] * B * phi[:, j]
        denom = max(norm(A * phi[:, j]), abs(lambda[j]) * norm(B * phi[:, j]), 1.0)
        @test norm(r[free]) / denom < tol
    end
end

@testset "solveEigen: single field, Dirichlet, symmetric, complex" begin
    line_mesh(n=8)
    P = Problem([Material("body")], type=:ScalarField, dim=1, fieldName=:p)
    nd = ndofs(P)
    A = spdiagm(0 => collect(1.0:nd))
    B = spdiagm(0 => ones(nd))
    KA = SystemMatrix(A, P)
    KB = SystemMatrix(B, P)
    bc = [BoundaryCondition("left"; p=0.0)]

    lambda, phi = solveEigen(KA, KB; support=bc, n=3, solver=:direct)
    @test length(lambda) == 3
    @test size(phi) == (nd, 3)
    fixed = constrainedDoFs(P, bc)
    free = freeDoFs(P, bc)
    @test all(iszero, phi[fixed, :])
    @test all(isapprox(norm(phi[:, j]), 1; rtol=1e-12) for j in 1:3)
    check_residual(A, B, lambda, phi, free)

    lambda_s, phi_s = solveEigen(Symmetric(KA), Symmetric(KB);
        support=bc, n=3, solver=:direct, normalize=:mass)
    @test sort(real.(lambda_s)) ≈ sort(real.(lambda))
    check_residual(A, B, lambda_s, phi_s, free)
    @test all(isapprox(real(dot(phi_s[:, j], B * phi_s[:, j])), 1; rtol=1e-12)
              for j in axes(phi_s, 2))

    # A real nonsymmetric matrix with imaginary eigenvalues ±i.
    D = spdiagm(0 => collect(5.0:nd+4.0))
    D[1,1] = 0.0
    D[2,2] = 0.0
    D[1,2] = -1.0
    D[2,1] = 1.0
    lambda_c, phi_c = solveEigen(SystemMatrix(D, P), KB;
        n=2, solver=:direct, shift=0.0, normalize=:unit)
    @test all(isapprox(abs(imag(z)), 1.0; rtol=1e-10) for z in lambda_c)
    @test eltype(phi_c) <: Complex
    check_residual(D, B, lambda_c, phi_c, collect(1:nd))
end

@testset "solveEigen: multifield and independent boundary conditions" begin
    line_mesh(n=6)
    P1 = Problem([Material("body")], type=:ScalarField, dim=1, fieldName=:p)
    P2 = Problem([Material("body")], type=:ScalarField, dim=1, fieldName=:T)
    nd = ndofs(P1)

    A11 = SystemMatrix(spdiagm(0 => collect(1.0:nd)), P1)
    A22 = SystemMatrix(spdiagm(0 => collect(11.0:10.0+nd)), P2)
    A12 = SystemMatrix(spzeros(nd, nd), P2, P1)
    A21 = SystemMatrix(spzeros(nd, nd), P1, P2)
    KA = SystemMatrix([A11 A12; A21 A22])

    B11 = SystemMatrix(spdiagm(0 => ones(nd)), P1)
    B22 = SystemMatrix(spdiagm(0 => ones(nd)), P2)
    B12 = SystemMatrix(spzeros(nd, nd), P2, P1)
    B21 = SystemMatrix(spzeros(nd, nd), P1, P2)
    KB = SystemMatrix([B11 B12; B21 B22])

    bc1 = BoundaryCondition("left"; problem=P1, p=0.0)
    bc2 = BoundaryCondition("left"; problem=P2, T=0.0)
    bc = [bc1, bc2]

    lambda, phi = solveEigen(KA, KB; support=bc, n=4, solver=:direct)
    @test length(lambda) == 4
    @test size(phi) == (2 * nd, 4)
    free, fixed, _ = LowLevelFEM.multifield_bc_data(KA, bc; nsteps=1)
    @test all(iszero, phi[fixed, :])
    check_residual(KA.A, KB.A, lambda, phi, free)
end

@testset "solveEigen: equality MPC" begin
    line_mesh(n=8)
    P = Problem([Material("body")], type=:ScalarField, dim=1, fieldName=:p)
    nd = ndofs(P)
    K = SystemMatrix(spdiagm(0 => collect(1.0:nd)), P)
    M = SystemMatrix(spdiagm(0 => ones(nd)), P)
    mpc = MPC(master="left", slave="right", field=P, p=true)
    lambda, phi = solveEigen(K, M; n=3, solver=:direct, mpc=[mpc])
    left = only(LowLevelFEM.nodesOnPhysicalGroup(P, "left"))
    right = only(LowLevelFEM.nodesOnPhysicalGroup(P, "right"))
    @test phi[left, :] ≈ phi[right, :]
    @test length(lambda) == 3
end

println("The tests passed if no @test failures were reported.")
