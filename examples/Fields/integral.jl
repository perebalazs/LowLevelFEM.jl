using LowLevelFEM
gmsh.initialize()
## 
gmsh.open("/home/perebal/Dokumentumok/GitHub/LowLevelFEM.jl/examples/Fields/cube1.geo")
##
mat = material("cube")
prob = Problem([mat])
##
load1 = load("cube", fx=1)
##
f1 = loadVector(prob, [load1])
display(f1)
##
gmsh.finalize()