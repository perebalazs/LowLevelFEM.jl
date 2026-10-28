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

## Transient integration

`HHT`, `CDM` and `FDM` share coordinate transformation, approximation and MPC
preparation. They return global nodal fields, with physical time samples
`t0 + (j-1)*Δt` for `j=1:n`. Input histories have either one constant column
or exactly `n` sampled columns; the solvers use column order and do not
interpolate the input time labels. `U0`/`X0` have one initial column or exactly
`n` columns; a complete history prescribes supported components while only its
first column initializes free components. `V0` supplies initial free velocities.

```julia
# Second-order dynamics: M*ü + C*u̇ + K*u = f
u, v = HHT(K, M, C, f, u0, v0, n, Δt;
           support=bc, mpc=ties, coordSys=cs, α=-0.1, δ=0.1)
u, v = CDM(K, M, C, f, u0, v0, n, Δt;
           support=bc, mpc=ties, coordSys=cs)

# First-order dynamics: C*ẋ + K*x = q
x = FDM(K, C, q, x0, n, Δt; support=bc, ϑ=0.5)
```

Omit the positional damping matrix in HHT/CDM for zero physical damping, or
supply it with `damping=C`. For the standard HHT parameter relation use
`δ=-α` (the historical independent `δ=0` default is retained). The historical positional boundary-list calls
remain available. Matrices and input fields are not modified.

HHT and FDM accept the same linear solver controls as `solveField`:
`solver`, `preconditioner`, `solveroptions` and `preconditioneroptions`.
Supported direct solvers are `:backslash`, `:lu`, `:lu_no_ordering`,
`:cholesky`, `:ldlt` and `:qr`. Iterative choices are `:cg` and `:gmres`.
The factorization or preconditioner of the effective free-coordinate matrix
is prepared once and reused across the entire time history. HHT uses
`M_ff + (1+α)*γ*Δt*C_ff + (1+α)*β*Δt^2*K_ff`; FDM uses
`C_ff + ϑ*Δt*K_ff`.

```julia
using IterativeSolvers, Preconditioners, IncompleteLU

# Declare all contributing matrices symmetric for CG + automatic ICH.
u, v = HHT(Symmetric(K), Symmetric(M), Symmetric(C), f, u0, v0, n, Δt;
           support=bc, solver=:cg, preconditioner=:auto,
           solveroptions=(reltol=1e-10, maxiter=1000),
           preconditioneroptions=(memory=3,))

# General/nonsymmetric systems: GMRES + ILU.
x = FDM(K, C, q, x0, n, Δt;
        support=bc, solver=:gmres, preconditioner=:ilu,
        solveroptions=(reltol=1e-10, maxiter=1000),
        preconditioneroptions=(τ=1e-3,))
```

`solver=:auto` is direct: ordinary matrices use Julia's factorization for
backslash, while explicitly symmetric systems try Cholesky with LU fallback.
Symmetry is preserved only when every contributing input matrix uses
`Symmetric(...)`; symmetric stiffness and mass do not imply symmetric damping.
CG additionally requires positive definiteness. `preconditioner=:auto` chooses
`:ich` for symmetric CG systems and `:ilu` for GMRES. `:none` and preconstructed
preconditioner objects are also accepted. Optional-package requirements and
native keyword validation match `solveField`.

For HHT with a nondiagonal mass matrix, the initial acceleration uses a separate
mass-system preparation with the same solver settings. This avoids an implicit
LU factorization when an iterative solver is requested. A diagonal mass still
uses direct division. CDM retains its diagonal/cached-mass solve.

Supports accept time callbacks as well as sampled scalar components:

```julia
bc = [BoundaryCondition("base"; field=U,
                        ux=(x,y,z,t)->0.001*sin(20t), uy=0.0)]
u, v = HHT(K, M, C, f, u0, v0, n, Δt; support=bc)
```

Displacement histories belong in `U0` (or `X0` for FDM), without additional
prescribed-history keywords:

- One column: initialize free DOFs; boundary values come from `support`.
- Exactly `n` columns: initialize free DOFs from column 1 and read the complete
  support history from the columns of `U0`/`X0`. Values in the BC objects are
  overridden. Future values on free DOFs are ignored and solved normally.
- Other column counts are rejected. `support` always selects constrained DOFs;
  a full initial field alone does not constrain any additional component.

```julia
# u_history is a nodal VectorField with n columns.
bc = [BoundaryCondition("base"; field=U, ux=0.0, uy=0.0)]
u, v = HHT(K, M, C, f, u_history, v0, n, Δt; support=bc)
x = FDM(K, C, q, x_history, n, Δt; support=bc)
```

Constrained velocity and acceleration follow the displacement history. HHT
estimates the initial derivatives from samples and then uses Newmark
kinematics; CDM uses finite differences. With only two samples the prescribed
history is linear and its initial acceleration is zero. FDM uses the complete
primary-variable increment, including capacity coupling. To specify a
velocity history, use `integrate(v_history)` (or `∫(v_history)`) to form `U0`
before calling the solver, choosing the required displacement integration
constant. `V0` controls only the initial free velocity.

For block systems all fields in a full history have the same column count.
For example, `U0=SystemVector([u_history, φ_history])` supplies a displacement
and rotation history together.

Ordinary loads, initial fields and prescribed histories use the local solver
basis. Wrap Cartesian data in `Global(...)`, just as for `solveField`.
Constrained DOFs use the full U0/X0 history when supplied, otherwise support
data. MPCs act on local components, including supports prescribed on slave nodes.

For block systems, use `SystemVector([field1,field2,...])` or a tuple in matrix
block order for loads and initial fields. Return values follow:

```julia
(u1,u2), (v1,v2) = HHT(Kblock, Mblock, Cblock, F, U0, V0, n, Δt;
                      support=bc, mpc=ties)
(u1,u2), (v1,v2) = CDM(Kblock, Mblock, Cblock, F, U0, V0, n, Δt;
                      support=bc, mpc=ties)
x1,x2 = FDM(Kblock, Cblock, F, X0, n, Δt; support=bc, mpc=ties)
```

Remote translation/rotation coupling uses the existing `rigidRotationMap`.
For `u_phys = u + R*φ`, form every physical operator with the same map:

```julia
R = rigidRotationMap(mpc_u, mpc_φ)
block(A) = SystemMatrix([A A*R; R'*A R'*A*R])
Kblock, Mblock, Cblock = block(Ku), block(Mu), block(Cu)
F = SystemVector([fu, R'*fu])

(u,φ), (v,ω) = HHT(Kblock, Mblock, Cblock, F, U0, V0, n, Δt;
                   support=bc, mpc=[mpc_u,mpc_φ])
u_phys = u + R*φ
v_phys = v + R*ω
```

The same block inputs work with CDM and FDM. Translation/rotation coupling in
mass and damping is essential: transforming stiffness alone leaves the
rotational inertia missing. No artificial independent mass on the auxiliary
rotation field is needed; rotational inertia follows from `R'*Mu*R`. Inactive
auxiliary DOFs are removed by the common MPC backend before integration.

Every block BC and MPC must explicitly identify its `field`. HHT and CDM
require a nonsingular mass matrix on the remaining free coordinates; they do
not solve arbitrary differential-algebraic systems. CDM uses diagonal division
when possible and a cached mass factorization otherwise. It does not silently
lump the mass after reduction, and its time step remains conditionally stable.

Reduced interpolation and SEM use the existing cached basis transformations.
Existing shared-backend limitations still apply: a coordinate basis cannot be
combined with reduced/spectral interpolation on the same field, and SEM plus
MPC on the same field is currently unsupported. Different block fields may
use these features separately.

`CDMaccuracyAnalysis` remains available. `HHTaccuracyAnalysis` includes a
physical damping ratio (`ξ`, a number or function of angular frequency) and
Rayleigh coefficients (`massDamping`, `stiffnessDamping`). For example:

```julia
x, radius = HHTaccuracyAnalysis(ωmin, ωmax, Δt, :SR;
                               α=-0.1, δ=0.1, ξ=0.02,
                               massDamping=0.0, stiffnessDamping=1e-4)
```

The scalar analysis describes modal/proportional damping. An arbitrary
nonproportional matrix is supported by the time integrator but cannot be
inferred from an undamped frequency interval. `FDMaccuracyAnalysis` returns
the spectral radius of the theta update, including for backward Euler.

```@docs
FDMaccuracyAnalysis
LowLevelFEM.prepare_transient_system
LowLevelFEM.reconstruct_transient_solution
```
