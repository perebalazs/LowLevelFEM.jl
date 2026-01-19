#include("../src/LowLevelFEM.jl")
using LowLevelFEM
using Test

@testset "Minimal cube problem" begin
    # Write your tests here.
    #@test LowLevelFEM.greet_your_package_name() == "Hello LowLevelFEM!"
    #@test LowLevelFEM.greet_your_package_name() != "Hello world!"
    #returntest = "test", 0, 0, 0
    #@test FEM.displacementConstraint("test", ux=0, uy=0, uz=0) == returntest
    #display(FEM.displacementConstraint("test", ux=0, uy=0, uz=0))
    try
        gmsh.initialize()
        if get(ENV, "CI", "false") == "true"
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.option.setNumber("General.Verbosity", 0)
        end
        gmsh.model.add("cube")
        b1 = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1, 1)
        pts = gmsh.model.getEntities(0)
        gmsh.model.occ.mesh.setSize(pts, 1.0)
        #gmsh.option.setNumber("Mesh.MeshSizeMin", 1.0)
        #gmsh.option.setNumber("Mesh.MeshSizeMax", 1.0)
        #gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(3)
        gmsh.model.addPhysicalGroup(3, [1], 10)
        gmsh.model.setPhysicalName(3, 10, "body")
        gmsh.model.addPhysicalGroup(2, [1], 11)
        gmsh.model.setPhysicalName(2, 11, "left")
        gmsh.model.addPhysicalGroup(2, [2], 12)
        gmsh.model.setPhysicalName(2, 12, "right")

        mat = material("body")
        prob = Problem([mat])

        supp = displacementConstraint("left", ux=0, uy=0, uz=0)
        load1 = load("right", fx=1)
        u = solveDisplacement(prob, load=[load1], support=[supp])
        S = solveStress(u)
        σ = probe(S, 0.5, 0.5, 0.5)
        @test all(!isnan, σ)
    finally
        gmsh.finalize()
    end
end
