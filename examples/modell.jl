import LowLevelFEM as FEM
using LowLevelFEM

gmsh.initialize()

gmsh.open("modell.geo")
mat = FEM.material("body", E=1.96e5, ν=0.25, ρ=7.874e-9)
problem = FEM.Problem([mat], type="PlaneStress")

supp = FEM.displacementConstraint("supp", ux=0, uy=0)
load = FEM.load("load", fx=-10)

K = FEM.stiffnessMatrix(problem)
M = FEM.massMatrix(problem, lumped=false)
f = FEM.loadVector(problem, [load])
C = K * 0

FEM.applyBoundaryConditions!(problem, K, M, C, f, [supp])

Tₘᵢₙ = FEM.smallestPeriodTime(K, M)

dof, dof = size(K)
u0 = zeros(dof)
v0 = zeros(dof)

display(1e-6/(Tₘᵢₙ / π))
u, v, t = FEM.CDM(K, M, C, f, u0, v0, 1e-6, Tₘᵢₙ / π)

uvec = FEM.showDoFResults(problem, u, t=t, "uvec", name="u(t)", visible=true)

gmsh.fltk.run()
gmsh.finalize()