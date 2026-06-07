using LowLevelFEM
gmsh.initialize()

gmsh.open("rectmesh.geo")
mat = material("body")
prob = Problem([mat], type=:PlaneStrain)

@time K = stiffnessMatrix(prob)
@time K = stiffnessMatrix(prob)
@time K = stiffnessMatrix(prob)

gmsh.finalize()
