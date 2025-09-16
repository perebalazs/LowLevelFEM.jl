using BenchmarkTools
using LowLevelFEM
gmsh.initialize()

cd("/home/perebal/Dokumentumok/GitHub/LowLevelFEM.jl/examples/multithread-test")
gmsh.open("cube.geo")

mat = material("body")
prob = Problem([mat])

@time K = LowLevelFEM.stiffnessMatrixSolidOld(prob)
@time K = LowLevelFEM.stiffnessMatrixSolid(prob)
@time K = LowLevelFEM.stiffnessMatrixSolidParallel(prob)

@btime K = LowLevelFEM.stiffnessMatrixSolidOld(prob)
@btime K = LowLevelFEM.stiffnessMatrixSolid(prob)
@btime K = LowLevelFEM.stiffnessMatrixSolidParallel(prob)

gmsh.finalize()
