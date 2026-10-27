# Assembly and Solvers

The modern workflow separates the weak form, assembled matrix/RHS, and solver.
All algebraic data remain accessible for inspection and custom formulations.

## A scalar diffusion problem

This example solves `-Δu=0` on the unit square with `u=0` on the left and
`u=1` on the right. The unconstrained top and bottom have zero natural flux.
The exact solution is `u(x,y)=x`.

```julia
# docs-check: scalar-diffusion
using LowLevelFEM

gmsh.initialize()
try
    structured_rect_mesh(n=2, order=1)
    mat = Material("body")
    P = Field([mat], type=:ScalarField, dim=2, fieldName=:u)
    K = ∫(Grad(P) ⋅ Grad(P))
    f = ∫(P ⋅ 0.0)
    bc = [BoundaryCondition("left", field=P, u=0.0),
          BoundaryCondition("right", field=P, u=1.0)]
    u = solveField(K, f; support=bc)
    exact = scalarField(P, "body", (x,y,z) -> x)
    @assert maximum(abs.(u.a - exact.a)) < 1e-10
finally
    gmsh.finalize()
end
```

In an existing session, manage Gmsh initialization/finalization once around the
whole workflow. A model and its fields must remain available while evaluating
or displaying results.

## Mixed problems

Assemble coupling blocks explicitly and preserve their field order:

```julia
K = SystemMatrix([Auu Aup; Aup' App])
F = SystemVector([fu, fp])
u, p = solveField(K, F; support=bc)
```

Here the four blocks/RHS fields are already assembled on compatible descriptors.
Bind each essential condition using `field=U` or `field=P`.
For lower-order pressure interpolation on a common higher-order mesh, define
its descriptor with `reducedOrder=true`; this differs from reduced quadrature.
A complete small coupled model is in the [matrix-level workflow](../explanations/matrix-level-workflow.md).

## Solver selection

The default general sparse solve uses backslash. `Symmetric(K)` allows an
attempted Cholesky factorization with LU fallback if the system is indefinite.
For a mixed saddle-point system, symmetry alone does not justify CG or Cholesky.

`solver=:cg`/`:gmres` require `using IterativeSolvers`. Automatic incomplete
Cholesky or ILU additionally requires `Preconditioners` or `IncompleteLU`.
Native iterative options are supplied through a `NamedTuple`, for example
`solveroptions=(reltol=1e-8,maxiter=1000)`. See [solveField](../reference/solution.md)
for all backends and coordinate/MPC restrictions.

For eigenproblems use `solveEigen(A,B; ...)`; interpretation and normalization
are described in [Solution API](../reference/solution.md).
For standard engineering problems, [linear mechanics](../reference/linear.md)
and [heat](../reference/heat.md) provide problem-oriented entry points.
These names remain while their backend migrates to DSL-backed wrappers;
see [API layers](../reference/legacy.md).
