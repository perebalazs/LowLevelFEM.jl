@testset "Plane stress validation" begin

    structured_rect_mesh()

    mat = Material("body")
    prob = Problem([mat], type=:PlaneStress)

    supp1 = displacementConstraint("left", ux=0)
    supp2 = displacementConstraint("bottom", uy=0)

    ld = load("right", fx=1)

    u = solveDisplacement(prob,
        load=[ld],
        support=[supp1, supp2]
    )

    u_FEM = probe(u, 1, 1, 0)

    ux_EX = 1 / mat.E
    uy_EX = -mat.ν * ux_EX

    @test isapprox(u_FEM[1], ux_EX; rtol=1e-2)
    @test isapprox(u_FEM[2], uy_EX; rtol=1e-2)

end