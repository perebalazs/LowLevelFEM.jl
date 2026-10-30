# Solution API

```@meta
CurrentModule = LowLevelFEM
```

## Matrices and right-hand sides

`SystemMatrix.A` exposes the sparse matrix. `SystemVector.a` and nodal field
`.a` expose RHS/solution matrices, with one column per independent RHS or mode.
For block systems, preserve the same `Problem` order in matrices and vectors.
`Symmetric(K, :U)` or `Symmetric(K, :L)` declares its authoritative triangle.
Solvers expand that triangle before coordinate, approximation or MPC projection,
so triangular storage is supported consistently across field, eigen and
transient solutions. The opposite triangle is intentionally ignored, following
Julia's `Symmetric` contract; projection roundoff follows the same convention.

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
depends on the chosen pencil. Dense `solver=:direct` has quadratic storage
cost; `dense_limit=1024^3` bytes guards an estimate of six full work arrays
before densification. Set a larger limit explicitly, or `nothing` to disable
the guard. This estimate is not a strict bound on backend peak memory.
For vibration, convert positive `λ=ω²` to Hz as
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

Normalized eigenvectors have their largest-magnitude component real and
positive; `normalize=:none` keeps the backend phase. Repeated eigenspaces can
still have different bases. The unshifted default `which=:SR` applies to
general pencils. For a symmetric positive spectrum, explicit `shift=0.0`
usually converges faster. `shift` is a strict numerical lower bound, not a
physical rigid-mode filter: choose a positive cutoff to exclude zero-mode
roundoff. Shift-invert regularization is placed below the requested bound and
requests additional candidates when necessary, retaining modes just above it.

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

## Transient integration

The [Time Integration guide](transient.md) contains the full HHT/CDM/FDM
workflow, keyword tables, support histories, solver choices and time-step
selection. The function references are [`HHT`](@ref), [`CDM`](@ref) and
[`FDM`](@ref).

| Method | System | Main choice |
|:--|:--|:--|
| HHT | Second order | Implicit dynamics and algorithmic damping. |
| CDM | Second order | Explicit displacement updates and a stability-limited step. |
| FDM | First order | Theta integration for diffusion/heat and other first-order systems. |

All three share coordinate, approximation and MPC preparation. Initial fields
have one column; `support` supplies each prescribed history once. Results are
global nodal fields at `t0+(j-1)*Δt`, including the initial sample.

## Operator modifications and nonlinear iteration

See [Operator updates and nonlinear iteration](@ref transient-callbacks) for
`modifyA2`, `modifyA1`, `modifyA0`, `modifyRHS`, callback state members and
nonlinear tolerances. Each callback defaults to `nothing`. HHT/FDM iterate
within a step; CDM keeps explicit displacement updates.

## First-order accuracy analysis

```@docs
FDMaccuracyAnalysis
```

## Transient preparation helpers

These lower-level helpers implement the shared preparation/reconstruction;
ordinary applications call HHT, CDM or FDM directly.

```@docs
LowLevelFEM.prepare_transient_system
LowLevelFEM.reconstruct_transient_solution
```
