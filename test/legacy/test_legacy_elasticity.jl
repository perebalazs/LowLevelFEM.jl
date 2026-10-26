@testset "Legacy plane-stress assembly and analytical solution" begin
    with_test_gmsh() do
        structured_rect_mesh(n=2, order=1)
        mat = Material("body")
        P = Problem([mat], type=:PlaneStress)
        K = stiffnessMatrix(P)
        @test size(K.A, 1) == size(K.A, 2)
        @test nnz(K.A) > 0
        @test norm(K.A - K.A') / norm(K.A) < 1e-10

        support = [displacementConstraint("left", ux=0),
                   displacementConstraint("bottom", uy=0)]
        u = solveDisplacement(P; load=[load("right", fx=1)], support=support)
        u_FEM = probe(u, 1, 1, 0)
        ux_EX = 1 / mat.E
        uy_EX = -mat.ν * ux_EX
        @test isapprox(u_FEM[1], ux_EX; rtol=1e-2)
        @test isapprox(u_FEM[2], uy_EX; rtol=1e-2)
    end
end

@testset "Legacy cube displacement and stress" begin
    with_test_gmsh() do
        structured_box_mesh(n=2, order=1)
        mat = material("body")
        P = Problem([mat])
        support = [displacementConstraint("left", ux=0, uy=0, uz=0)]
        u = solveDisplacement(P; load=[load("right", fx=1)], support=support)
        stress = solveStress(u)
        @test all(isfinite, probe(stress, 0.5, 0.5, 0.5))
    end
end
