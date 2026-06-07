using LinearAlgebra, SparseArrays
import LowLevelFEM as FEM
using LowLevelFEM

gmsh.initialize()

gmsh.open("foam2D.geo")
problem = FEM.Problem(type="PlaneStress")

#gmsh.fltk.run()

supp = FEM.displacementConstraint("supp", ux=0, uy=0)
load = FEM.load("load", fx=-1)

K = FEM.stiffnessMatrix(problem)
M = FEM.massMatrix(problem)
C = K * 0
dropzeros!(C)
f = FEM.loadVector(problem, [load])
dof = problem.non * problem.dim
u0 = zeros(dof)
v0 = zeros(dof)

FEM.applyBoundaryConditions!(problem, K, M, C, f, [supp])

Tₘᵢₙ = FEM.smallestPeriodTime(K, M)

u, v, t = FEM.CDM(K, M, C, f, u0, v0, Tₘᵢₙ * 3, Tₘᵢₙ / π)

f *= 0
nst = size(u, 2)
u, v, t = FEM.CDM(K, M, C, f, u[:,nst], v[:,nst], Tₘᵢₙ * 500, Tₘᵢₙ / π)


#q = FEM.solveDisplacement(K, f)
S = FEM.solveStress(problem, u)

#u = FEM.showDoFResults(problem, u, "uvec", t=t, name="uvec", visible=false)
v = FEM.showDoFResults(problem, v, "vvec", t=t, name="vvec", visible=false)
#ux = FEM.showDoFResults(problem, q, "ux", name="ux", visible=false)
#uy = FEM.showDoFResults(problem, q, "uy", name="uy", visible=false)
#uz = FEM.showDoFResults(problem, q, "uz", name="uz", visible=false)

#s = FEM.showStressResults(problem, S, "s", t=t, name="σ", visible=true, smooth=true)
s = FEM.showStressResults(problem, S, "sx", t=t, name="σ", visible=true, smooth=true)
#sx = FEM.showStressResults(problem, S, "sx", name="σx", visible=false, smooth=true)
#sy = FEM.showStressResults(problem, S, "sy", name="σy", visible=false, smooth=true)
#sz = FEM.showStressResults(problem, S, "sz", name="σz", visible=false, smooth=true)
#sxy = FEM.showStressResults(problem, S, "sxy", name="τxy", visible=false, smooth=true)
#syz = FEM.showStressResults(problem, S, "syz", name="τyz", visible=false, smooth=true)
#szx = FEM.showStressResults(problem, S, "szx", name="τzx", visible=false, smooth=true)

#FEM.plotOnPath(problem, "path", sx, 100, name="σx", visible=false);
#FEM.plotOnPath(problem, "path", sxy, 100, name="τxy", visible=false);
#FEM.plotOnPath(problem, "path", ux, 100, name="ux", visible=false);

gmsh.fltk.run()
gmsh.finalize()
