import LowLevelFEM as FEM
using LowLevelFEM
using SparseArrays
using Plots

gmsh.initialize()
gmsh.open("DINAMIC_BIG_2.geo");
#mat = FEM.material("body")
mat1 = FEM.material("Fe", E=1.96e5, ν=0.25, ρ=7.874e-9)
mat2 = FEM.material("Al", E=7e4, ν=0.334, ρ=2.7e-9)
problem = FEM.Problem([mat1, mat2], type="PlaneStress")
#problem = FEM.Problem([mat], type="PlaneStress")
load = FEM.load("right", fx=-1);
print("K..")
K = FEM.stiffnessMatrix(problem);
println("OK")
print("M..")
M = FEM.massMatrix(problem, lumped=true);
println("OK")
print("f..")
f = FEM.loadVector(problem, [load]);
println("OK")
dof = problem.non * problem.dim;

u0 = zeros(dof);
v0 = zeros(dof);

C = spzeros(dof, dof);
C = 1e-8 * K;
print("Tₘᵢₙ..")
Tₘᵢₙ = FEM.smallestPeriodTime(K, M)
println("OK")
n = 3000; # egy lépéssel kevesebb lesz végül az átfedés miatt!

u = zeros(dof, n - 1);
v = zeros(dof, n - 1);
t = zeros(n - 1);

terhelt = 4;
terheletlen = n - terhelt;

print("u_1, v_1, t_1..")
u_1, v_1, t_1 = FEM.CDM(K, M, C, f, u0, v0, terhelt * Tₘᵢₙ / π / 1.3, Tₘᵢₙ / π / 1.3);
println("OK")

u_1_kezdo = u_1[:, terhelt];
v_1_kezdo = v_1[:, terhelt];

#display(K)
#display(M)
#display(C)
#display(f)
#display(u_1_kezdo)
#display(v_1_kezdo)
#display(terheletlen)
#display(Tₘᵢₙ)

print("u_2, v_2, t_2..")
u_2, v_2, t_2 = FEM.CDM(K, M, C, f * 0, u_1_kezdo, v_1_kezdo, terheletlen * Tₘᵢₙ / π / 1.3, Tₘᵢₙ / π / 1.3);
println("OK")



t_eltolas = t_1[terhelt];
t_2 = t_2 .+ t_eltolas;

u[:, 1:terhelt] = u_1;
u[:, terhelt:n-1] = u_2;

v[:, 1:terhelt] = v_1;
v[:, terhelt:n-1] = v_2;

t[1:terhelt] = t_1;
t[terhelt:n-1] = t_2;
print("S..")
S = FEM.solveStress(problem, u);
println("OK")
print("s..")
s = FEM.showStressResults(problem, S, "s", t=t, name="σred", visible=true, smooth=false)
println("OK")
print("uvec..")
uvec = FEM.showDoFResults(problem, u, "uvec", t=t, name="u", visible=false)
println("OK")
print("vvec..")
vvec = FEM.showDoFResults(problem, v, "vvec", t=t, name="v", visible=false)
println("OK")
#Wi = [(u[:, i] - u[:, i-1])' * (i > 5 ? f * 0 : f) for i in 2:n-1]
Wi = [(u[:, i] - u[:, i-1])' * f for i in 2:4]
W = sum(Wi)
U = u[:, n-1]' * K * u[:, n-1] / 2 + v[:, n-1]' * M * v[:, n-1] / 2
nn = 4    # vagy 5
U = u[:, nn]' * K * u[:, nn] / 2 + v[:, nn]' * M * v[:, nn] / 2
h = 2.5
np = 100
q = 0
X = 7.4999
#display(t)
print("q_1..")
q_1 = zeros(n - 1, 1);
for j ∈ 2:terhelt+2
    for i ∈ 1:np
        ΔA = h / np
        y = (h / np) * (i - 1) + (h / np / 2)
        Δt = t[j] - t[j-1]
        ss = reshape(gmsh.view.probe(s, X, y, 0, j - 1)[1], 3, 3)
        vv = gmsh.view.probe(vvec, X, y, 0, j - 1)[1]
        q += [1, 0, 0]' * ss * vv * ΔA * Δt
    end
    q_1[j, 1] = q
end
#display(q_1)
println("OK")
q
plot(t, q_1)

# Kompenzálás csillapításból
C_E = 0;
for i ∈ 1:terhelt
    C_E += v[:, i]' * C * v[:, i]
end
#display(U)
#display(W)
#display(q)
(U - W) / U * 100
(U - q) / U * 100
(W - q) / W * 100
((U + C_E) - W) / U * 100
(U + C_E - q) / U * 100

h = 2.5
np = 100
q = 0
X = 5
#display(t)
q_2 = zeros(n - 1, 1);
print("q_2..")
for j ∈ 2:n-1
    for i ∈ 1:np
        ΔA = h / np
        y = (h / np) * (i - 1) + (h / np / 2)
        Δt = t[j] - t[j-1]
        ss = reshape(gmsh.view.probe(s, X, y, 0, j - 1)[1], 3, 3)
        vv = gmsh.view.probe(vvec, X, y, 0, j - 1)[1]
        q += [1, 0, 0]' * ss * vv * ΔA * Δt
    end
    q_2[j, 1] = q
end
println("OK")
#display(q_2)
q
plot!(t, q_2)
# kérdéses, hogy miért nem csökken vissza, amikor a 2. fronton túlhalad a sebesség hullám.
h = 2.5
np = 100
q = 0
X = 0.001
#display(t)
q_3 = zeros(n - 1, 1);
print("q_3..")
for j ∈ 2:n-1
    for i ∈ 1:np
        ΔA = h / np
        y = (h / np) * (i - 1) + (h / np / 2)
        Δt = t[j] - t[j-1]
        ss = reshape(gmsh.view.probe(s, X, y, 0, j - 1)[1], 3, 3)
        vv = gmsh.view.probe(vvec, X, y, 0, j - 1)[1]
        q += [1, 0, 0]' * ss * vv * ΔA * Δt
    end
    q_3[j, 1] = q
end
println("OK")
#display(q_3)
q
plot!(t, q_3, label=["first surface" "second surface" "third surface"], xlabel="t [s]", ylabel="Work [mJ]")
print("pdf..")
Plots.pdf("surfaces.pdf")
println("OK")
gmsh.fltk.run()
gmsh.finalize()