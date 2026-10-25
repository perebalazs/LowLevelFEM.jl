using Test
using SparseArrays
using LinearAlgebra
using LowLevelFEM

# Tests the algebraic solver independently of Gmsh / problem assembly.
# Run with: include("eigen_api_tests.jl")

const solve_reduced = LowLevelFEM.solve_eigen_reduced

function diagonal_pencil(values)
    A = spdiagm(0 => Float64.(values))
    B = spdiagm(0 => ones(Float64, length(values)))
    return (; A, B, symmetricA=true, symmetricB=true, uploA=:U, uploB=:U)
end

function check_eigenpairs(P, λ, ϕ; tol=1e-8)
    @test length(λ) == size(ϕ, 2)
    @test issorted(real.(λ))
    @test all(norm(P.A * ϕ[:, i] - λ[i] * P.B * ϕ[:, i]) /
              (norm(P.A * ϕ[:, i]) + abs(λ[i]) * norm(P.B * ϕ[:, i])) < tol
              for i in eachindex(λ))
end

@testset "solveEigen: ordered generalized eigenpairs" begin
    P = diagonal_pencil(1:40)

    @testset "Direct: lowest, threshold, native sigma" begin
        λ, ϕ = solve_reduced(P; n=4, solver=:direct)
        @test λ ≈ [1, 2, 3, 4]
        check_eigenpairs(P, λ, ϕ)

        λ, ϕ = solve_reduced(P; n=4, shift=12.0, solver=:direct)
        @test λ ≈ [13, 14, 15, 16]
        check_eigenpairs(P, λ, ϕ)

        λ, ϕ = solve_reduced(P; n=3, sigma=12.2, which=:LM, solver=:direct)
        @test λ ≈ [11, 12, 13]
        check_eigenpairs(P, λ, ϕ)

        λ, _ = solve_reduced(P; n=3, which=:LM, solver=:direct)
        @test λ ≈ [38, 39, 40]
    end

    @testset "ARPACK: one-sided shift and native options" begin
        λ, ϕ = solve_reduced(P; n=4, shift=12.0, solver=:arpack,
                             solveroptions=(;ncv=24, tol=1e-11))
        @test λ ≈ [13, 14, 15, 16] atol=1e-6
        check_eigenpairs(P, λ, ϕ)

        λ, ϕ = solve_reduced(P; n=3, sigma=12.2, which=:LM,
                             solver=:arpack, solveroptions=(;ncv=24))
        @test λ ≈ [11, 12, 13] atol=1e-6
        check_eigenpairs(P, λ, ϕ)
    end

    @testset "Validation" begin
        @test_throws ArgumentError solve_reduced(P; n=4, shift=12.0, sigma=12.0)
        @test_throws ArgumentError solve_reduced(P; n=4, shift=12.0, which=:LM)
        @test_throws ArgumentError solve_reduced(P; n=4, which=:UNKNOWN)
        @test_throws ArgumentError solve_reduced(P; n=4, solver=:unknown)
        @test_throws ArgumentError solve_reduced(P; n=4, solveroptions=(;invalid=4))
        @test_throws ArgumentError solve_reduced(P; n=4, solveroptions=(;ritzvec=false))
    end

    @testset "Complex direct: ordered by real part" begin
        A = spdiagm(0 => ComplexF64[3+2im, -2+1im, 5-3im, 1-1im])
        B = spdiagm(0 => ones(ComplexF64, 4))
        Q = (;A, B, symmetricA=false, symmetricB=false, uploA=:U, uploB=:U)
        λ, ϕ = solve_reduced(Q; n=2, shift=0.0, solver=:direct)
        @test real.(λ) ≈ [1, 3]
        check_eigenpairs(Q, λ, ϕ)
        @test_throws ArgumentError solve_reduced(Q; n=2, shift=0.0, solver=:arpack)
    end
end
