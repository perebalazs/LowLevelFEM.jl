using BenchmarkTools
using LowLevelFEM
gmsh.initialize()

cd("/home/perebal/Dokumentumok/GitHub/LowLevelFEM.jl/examples/multithread-test")
gmsh.open("cube.geo")

mat = material("body")
prob = Problem([mat])

@time K = stiffnessMatrix(prob)

@btime K = stiffnessMatrix(prob)
@btime K = stiffnessMatrix(prob, forceOneThread=true)

gmsh.finalize()
