using LinearAlgebra, SparseArrays
include("./LowLevelFEM.jl")
import .LowLevelFEM as FEM
import gmsh_jll
include(gmsh_jll.gmsh_api)
gmsh.initialize()

gmsh.open("foam2D.geo")
problem = FEM.Problem(type="PlaneStress")

#gmsh.fltk.run()

supp = FEM.displacementConstraint("supp", ux=0, uy=0)
load = FEM.load("load", fx=1)

K = FEM.stiffnessMatrix(problem)
display(K)
f = FEM.loadVector(problem, [load])
display(size(K))

FEM.applyBoundaryConditions!(problem, K, f, [supp])

q = FEM.solveDisplacement(K, f)
S = FEM.solveStress(problem, q)

u = FEM.showDoFResults(problem, q, "uvec", name="uvec", visible=false)
ux = FEM.showDoFResults(problem, q, "ux", name="ux", visible=false)
uy = FEM.showDoFResults(problem, q, "uy", name="uy", visible=false)
uz = FEM.showDoFResults(problem, q, "uz", name="uz", visible=false)

s = FEM.showStressResults(problem, S, "s", name="σ", visible=true, smooth=true)
sx = FEM.showStressResults(problem, S, "sx", name="σx", visible=false, smooth=true)
sy = FEM.showStressResults(problem, S, "sy", name="σy", visible=false, smooth=true)
sz = FEM.showStressResults(problem, S, "sz", name="σz", visible=false, smooth=true)
sxy = FEM.showStressResults(problem, S, "sxy", name="τxy", visible=false, smooth=true)
syz = FEM.showStressResults(problem, S, "syz", name="τyz", visible=false, smooth=true)
szx = FEM.showStressResults(problem, S, "szx", name="τzx", visible=false, smooth=true)

#FEM.plotOnPath(problem, "path", sx, 100, name="σx", visible=false);
#FEM.plotOnPath(problem, "path", sxy, 100, name="τxy", visible=false);
#FEM.plotOnPath(problem, "path", ux, 100, name="ux", visible=false);

gmsh.fltk.run()
gmsh.finalize()
