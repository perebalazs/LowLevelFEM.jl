using LowLevelFEM

structured_rect_mesh()

mat = Material("body")
prob = Problem([mat], type=:PlaneStress)

supp1 = displacementConstraint("left", ux=0)
supp2 = displacementConstraint("bottom", uy=0)

ld = load("right", fx=1)

u = solveDisplacement(prob, load=[ld], support=[supp1,supp2])

E = solveStrain(u)
S = solveStress(u)

showDoFResults(u, name="u")
showStrainResults(E, :y)
showStressResults(S, :xy)

u_FEM = probe(u, 1,1,0)
u_FEM[1]
u_FEM[2]

ux_EX = 1 * 1 / mat.E / 1
uy_EX = -mat.ν * ux_EX

if u_FEM[1] ≈ ux_EX
    println("PlaneStress: ux -> OK")
end

if u_FEM[2] ≈ uy_EX
    println("PlaneStress: uy -> OK")
end