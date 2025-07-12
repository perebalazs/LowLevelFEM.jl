cd(@__DIR__)

import LowLevelFEM as FEM
using LowLevelFEM

gmsh.initialize()

gmsh.open("henger.geo")

mat = FEM.material("body", E=2.e5, ν=0.3)
problem = FEM.Problem([mat], type=:Solid)

supp = FEM.displacementConstraint("outer", ux=0, uy=0, uz=0)

p=10
f1(x,y,z) = p*(x/sqrt(x^2+y^2))
f2(x,y,z) = -p*(y/sqrt(x^2+y^2))

load = FEM.load("inner", fx=f2, fy=f1)

q = FEM.solveDisplacement(problem, [load], [supp])
S = FEM.solveStress(problem, q)

u = FEM.showDoFResults(problem, q, :uvec)
ux = FEM.showDoFResults(problem, q, :ux)
uy = FEM.showDoFResults(problem, q, :uy)

s = FEM.showStressResults(problem, S, :s, visible=true, smooth=true)
sx = FEM.showStressResults(problem, S, :sx, name="σx", visible=false, smooth=true)
sy = FEM.showStressResults(problem, S, :sy, name="σy", visible=false, smooth=true)
sxy = FEM.showStressResults(problem, S, :sxy, name="τxy", visible=false, smooth=true)

gmsh.fltk.run()
gmsh.finalize()