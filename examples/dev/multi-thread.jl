using LowLevelFEM, LinearAlgebra

println("Szálak száma: ", Threads.nthreads())
println("LinAlg szálak száma: ", LinearAlgebra.BLAS.get_num_threads())

structured_box_mesh(n=10, order=2)

mat = Material("body")
Pu = Problem([mat], type=:VectorField, dim=3, field=:u)

prob = Problem([mat])
println("Legacy stiffness matrix:")
stiffnessMatrix(prob)
GC.gc()
@time K1 = stiffnessMatrix(prob)

μ = mat.μ
λ = mat.λ
D = [λ+2μ λ λ 0 0 0; λ λ+2μ λ 0 0 0; λ λ λ+2μ 0 0 0; 0 0 0 μ 0 0; 0 0 0 0 μ 0; 0 0 0 0 0 μ]

println("Stiffness matrix with DSL:")
∫(SymGrad(Pu) ⋅ D ⋅ SymGrad(Pu))
GC.gc()
@time K2 = ∫(SymGrad(Pu) ⋅ D ⋅ SymGrad(Pu))

print("Error: ")
err = norm(K1.A - K2.A) / norm(K1.A)
println(100err, "%")
println()

structured_rect_mesh(x0=10.0, n=50, order=2)

prob = Problem([mat], type=:AxiSymmetric)

println("Legacy axisymmetrix problem:")
stiffnessMatrix(prob)
GC.gc()
@time K1 = stiffnessMatrix(prob)

Pu = Problem([mat], type=:VectorField, dim=2, field=:u)

E = mat.E
ν = mat.ν

r = ScalarField(Pu, "body", (x, y, z)->x)
A1 = [1 0 0; 0 0 0; 0 1 0; 0 0 1]
A2 = [0 0; 1/r 0; 0 0; 0 0]
B = A1 ⋅ SymGrad(Pu) + A2 ⋅ Pu
D = E / (1+ν) / (1-2ν) * [1-ν ν ν 0; ν 1-ν ν 0; ν ν 1-ν 0; 0 0 0 (1-2ν)/2]

println("Axisymmetric problem with compound DSL:")
∫(B' ⋅ D ⋅ B * (2π*r))
GC.gc()
@time K2 = ∫(B' ⋅ D ⋅ B * (2π*r))

print("Error: ")
err = norm(K1.A - K2.A) / norm(K1.A)
println(100err, "%")
println()
println("LinAlg szálak száma: ", LinearAlgebra.BLAS.get_num_threads())