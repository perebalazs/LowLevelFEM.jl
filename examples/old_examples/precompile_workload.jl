using LowLevelFEM
using gmsh

# kis segédfüggvény: létrehoz egy triviális problem objektumot
function make_test_problem(; dim::Int, type::Symbol, phname::String="Mat1")
    gmsh.initialize()
    gmsh.model.add("test")

   dim = 0
    if type == :Solid
         dim = 3
        gmsh.model.occ.addRectangle(0, 0, 0, 1, 1, 1)
        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(2, [1], 1)
        gmsh.model.setPhysicalName(2, 1, phname)
        gmsh.model.mesh.generate(2)
    elseif type == :PlaneStress || type == :AxiSymmetric
         dim = 2
        gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1, 1)
        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(3, [1], 1)
        gmsh.model.setPhysicalName(3, 1, phname)
        gmsh.model.mesh.generate(3)
    else
        error("type=$type nem támogatott")
    end

    # egyszerű homogén anyag
    mat = (phName=phname, E=210e9, ν=0.3, ρ=7800.0)

    Problem([mat], type=type)
    
end

# tipikus problémák
prob2d_ps  = make_test_problem(dim=2, type=:PlaneStress)
#prob2d_axi = make_test_problem(dim=2, type=:AxiSymmetric)
#prob3d     = make_test_problem(dim=3, type=:Solid)

# előfordítandó hívások
stiffnessMatrixSolid(prob2d_ps)
#stiffnessMatrixAxi(prob2d_axi)
#stiffnessMatrixSolid(prob3d)

massMatrixSolid(prob2d_ps)
massMatrixSolid(prob2d_axi)
massMatrixSolid(prob3d)

gmsh.finalize()

