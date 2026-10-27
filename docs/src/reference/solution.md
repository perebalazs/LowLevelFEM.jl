# Solution API

```@meta
CurrentModule = LowLevelFEM
```

## Matrices and right-hand sides

`SystemMatrix.A` exposes the sparse matrix. `SystemVector.a` and nodal field
`.a` expose RHS/solution matrices, with one column per independent RHS or mode.
For block systems, preserve the same `Problem` order in matrices and vectors.

```@docs
SystemMatrix
SystemVector
subSystemMatrix
```

## Linear field solution

```@docs
solveField
Global
```

The tutorial [matrix-level workflow](../explanations/matrix-level-workflow.md)
includes single-field, coupled-field and eigen examples. Direct solvers are
available in the base package; iterative solvers and preconditioners require
explicitly loaded optional packages described by `solveField`.

## Eigenproblems

`solveEigen(A,B)` solves `A*ϕ = λ*B*ϕ`; the physical interpretation of `λ`
depends on the chosen pencil. For vibration, convert positive `λ=ω²` to Hz as
`sqrt(λ)/(2π)`. A buckling multiplier is not a frequency.

```@docs
solveEigen
solveEigenProblem
eigenValue
eigenVector
showEigenResults
getEigenValues
LowLevelFEM.mapEigenValues
```

`mapEigenValues` currently accepts `VectorField` eigenresults; it does not have
scalar/tensor overloads. For other field types, extract values with `eigenValue`
or request `output=:raw`.

See [Linear Buckling and Prestressed Natural Frequencies](../tutorials/buckling-and-prestressed-modes.md)
for a complete mechanical application. Engineering mechanical eigen helpers and their existing result containers
are documented under [Linear API](linear.md).

## Interpolation-space transformations

Field storage remains in Gmsh/Lagrange coordinates even when the solver uses a
reduced or spectral basis. These are advanced transformation APIs; the ordinary
solver constructs and applies its required transformations automatically.

```@docs
reductionMatrices
reduced_bc_data
spectralTransformation
```

Reduced quadrature (`reduced(expr)`) is distinct from reduced interpolation
(`reducedOrder=true`). Spectral support currently covers lines and complete
quadrangles; consult `spectralTransformation` for unsupported combinations.
