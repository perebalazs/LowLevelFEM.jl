# Run within your LowLevelFEM development environment after replacing the sources.
# All docstrings and test comments are in English.
using Test, LinearAlgebra, SparseArrays, LowLevelFEM

@testset "Field-valued solveEigen API" begin
    structured_rect_mesh(n=2)
    mat = Material("body")
    P = Problem([mat], type=:VectorField, dim=2, fieldName=:u)
    nd = ndofs(P)
    K = SystemMatrix(spdiagm(0 => collect(1.0:nd)), P)
    M = SystemMatrix(spdiagm(0 => ones(nd)), P)

    u = solveEigen(K, M; nev=4, solver=:direct)
    @test u isa VectorField
    @test u.kind === :eigen_real
    @test u.type === :v2D
    @test u.nsteps == 4
    @test isNodal(u)
    @test size(u.a) == (nd, 4)
    @test u.t ≈ collect(1.0:4.0)
    @test eigenValue(u; n=3) ≈ 3.0
    u3 = eigenVector(u; n=3)
    @test u3 isa VectorField
    @test u3.kind === :eigen_real
    @test u3.nsteps == 1
    @test isNodal(u3)
    @test u3.a == u.a[:, 3:3]
    @test eigenValue(u3) ≈ 3.0
    @test copy(u).kind === :eigen_real
    @test nodesToElements(u3).kind === :eigen_real
    @test elementsToNodes(nodesToElements(u3)).kind === :eigen_real

    ui = solveEigen(K, M; nev=4, solver=:direct, part=:imag)
    @test ui.kind === :eigen_imag
    @test all(iszero, ui.t)
    ur, ui = solveEigen(K, M; nev=4, solver=:direct, part=:both)
    @test ur.kind === :eigen_real && ui.kind === :eigen_imag
    @test ur.t == u.t

    λ, ϕ = solveEigen(K, M; nev=4, solver=:direct,
                      output=:raw, part=:invalid)
    @test λ ≈ u.t
    @test size(ϕ) == size(u.a)

    normal = VectorField(Matrix{Float64}[], zeros(nd, 1), [0.0],
                         Int[], 1, :v2D, P)
    @test normal.kind === :field
    @test_throws ArgumentError eigenValue(normal)
    @test_throws ArgumentError eigenVector(normal)
    @test_throws ArgumentError solveEigen(K, M; output=:unknown)
    @test_throws ArgumentError solveEigen(K, M; part=:unknown)

    # These use the ordinary Gmsh field output and should create View tags.
    @test showEigenResults(u; visible=false) isa Integer
    @test showEigenResults(u3; visible=false) isa Integer
    @test showModalResults(u; visible=false) isa Integer
    @test showBucklingResults(u; visible=false) isa Integer

    # A skew 2x2 block creates genuinely complex eigenpairs from real matrices.
    A = spdiagm(0 => collect(2.0:nd+1.0))
    A[1, 1] = 0.0; A[2, 2] = 0.0
    A[1, 2] = -1.0; A[2, 1] = 1.0
    C = SystemMatrix(A, P)
    zr, zi = solveEigen(C, M; nev=2, solver=:direct, part=:both)
    @test zr isa VectorField && zi isa VectorField
    @test any(x -> abs(x) > 0.5, zi.t)
    @test any(x -> abs(x) > 0.5, zi.a)

    # Mixed scalar-vector multifield matrices, no physical coupling required.
    Ps = Problem([mat], type=:ScalarField, dim=2, fieldName=:s)
    Pv = Problem([mat], type=:VectorField, dim=2, fieldName=:v)
    ns, nv = ndofs(Ps), ndofs(Pv)
    Ks = SystemMatrix(spdiagm(0 => collect(1.0:ns)), Ps)
    Kv = SystemMatrix(spdiagm(0 => collect(100.0:(99+nv))), Pv)
    Ms = SystemMatrix(spdiagm(0 => ones(ns)), Ps)
    Mv = SystemMatrix(spdiagm(0 => ones(nv)), Pv)
    Zsv = SystemMatrix(spzeros(ns, nv), Pv, Ps)
    Zvs = SystemMatrix(spzeros(nv, ns), Ps, Pv)
    Kb = SystemMatrix([Ks Zsv; Zvs Kv])
    Mb = SystemMatrix([Ms Zsv; Zvs Mv])

    s, v = solveEigen(Kb, Mb; nev=4, solver=:direct)
    @test s isa ScalarField && v isa VectorField
    @test s.kind == :eigen_real && v.kind == :eigen_real
    @test s.t == v.t
    @test size(s.a) == (ns, 4)
    @test size(v.a) == (nv, 4)
    sr, si, vr, vi = solveEigen(Kb, Mb; nev=4, solver=:direct, part=:both)
    @test sr isa ScalarField && si isa ScalarField
    @test vr isa VectorField && vi isa VectorField
    @test eigenValue(sr; n=2) ≈ eigenValue(vr; n=2)

    # Preserve old raw access for internal solvers.
    oldλ, oldϕ = solveEigenProblem(K, M; n=3, directSolver=true)
    @test length(oldλ) == 3
    @test size(oldϕ) == (nd, 3)
end
