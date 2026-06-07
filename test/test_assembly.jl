@testset "Assembly" begin
    gmsh.initialize()
    try
        structured_rect_mesh()

        mat = Material("body")
        P = Problem([mat], type=:PlaneStress)

        K = stiffnessMatrix(P)

        @test size(K.A,1) == size(K.A,2)
        @test nnz(K.A) > 0
        @test norm(K.A - K.A') / norm(K.A) < 1e-10

    finally
        gmsh.finalize()
    end
end