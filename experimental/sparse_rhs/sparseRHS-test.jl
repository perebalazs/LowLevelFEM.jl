using SparseArrays, LinearAlgebra
using .SparseSolveAuto

A = sprand(1000, 1000, 0.001) + I
b = sprand(1000, 0.001)
C = sprand(1000, 10, 0.001)

F = lu(A)

x = SparseSolveAuto.solve_sparse_rhs(F, b)
println("x típusa: ", typeof(x))

X = SparseSolveAuto.solve_sparse_rhs(F, C)
println("X típusa: ", typeof(X))

