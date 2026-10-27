# Modal and Buckling Analysis

Use [Linear Buckling and Prestressed Natural Frequencies](buckling-and-prestressed-modes.md)
for the complete DSL tutorial and executable notebook. It covers static preload,
geometric stiffness, eigenvalue mappings and interpretation of mode shapes.

The generic API is [solveEigen](../reference/solution.md). Engineering mechanical
helpers are documented under [Linear mechanics](../reference/linear.md).
This page retains the previous tutorial URL.

## A small exact discrete eigenproblem

The following base-package example needs no ARPACK. Two linear elements on a
unit interval leave one unconstrained interior DOF. Its stiffness is `4` and
mass is `1/3`, giving the discrete eigenvalue `12`.

```julia
# docs-check: scalar-eigen
using LowLevelFEM, LinearAlgebra

gmsh.initialize()
try
    line_mesh(n=2, order=1)
    P = Field([Material("body")], type=:ScalarField, dim=1, fieldName=:p)
    K = ∫(Grad(P) ⋅ Grad(P))
    M = ∫(P ⋅ P)
    bc = [BoundaryCondition("left", field=P, p=0.0),
          BoundaryCondition("right", field=P, p=0.0)]
    modes = solveEigen(K, M; n=1, solver=:direct, normalize=:mass, support=bc)
    @assert isapprox(eigenValue(modes), 12.0; atol=1e-10)
    phi = eigenVector(modes)
    @assert isapprox(only(phi.a' * M.A * phi.a), 1.0; atol=1e-10)
finally
    gmsh.finalize()
end
```

The continuum first eigenvalue is `π²`; refinement approaches it. Passing the
assertions checks this discretization and normalization, not mesh convergence.
