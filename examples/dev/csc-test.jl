using LowLevelFEM, LinearAlgebra, SparseArrays
structured_rect_mesh(x0=10.0, n=200, order=2)
mat = Material("body")
Pu = Problem([mat], type=:VectorField, dim=2, field=:u)

E = mat.E
ν = mat.ν

r = ScalarField(Pu, "body", (x, y, z)->x)
A1 = [1 0 0; 0 0 0; 0 1 0; 0 0 1]
A2 = [0 0; 1/r 0; 0 0; 0 0]
B = A1 ⋅ SymGrad(Pu) + A2 ⋅ Pu
D = E / (1+ν) / (1-2ν) * [1-ν ν ν 0; ν 1-ν ν 0; ν ν 1-ν 0; 0 0 0 (1-2ν)/2]


expr = B' ⋅ D ⋅ B

K = ∫(
    expr * (2π*r);
    Ω="body",
    assembly=:ijv,
    multithread=false
)

GC.gc()
@time K = ∫(
    expr * (2π*r);
    Ω="body",
    assembly=:csc_search,
    multithread=false
)

println("DOF: ", size(K.A))
println("nnz: ", nnz(K.A))
