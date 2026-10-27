# Engineering API, DSL and Legacy Assembly

```@meta
CurrentModule = LowLevelFEM
```

LowLevelFEM provides two public ways to formulate a model:

- **Engineering API:** problem-oriented names such as `stiffnessMatrix`,
  `massMatrix`, `loadVector`, `solveDisplacement` and `solveTemperature`.
  These names are intended to remain, with broadly the same argument patterns.
- **Weak-form DSL:** explicit operators, coefficients, matrix/load assembly and
  generic field solvers for custom formulations and multiphysics.

The intended architecture uses the DSL engine behind the engineering functions:

```text
engineering API (materials, loads, supports, problem-specific helpers)
                            ↓ wrappers
             weak-form DSL / assembly engine
                            ↓
                  matrices, solvers, fields
```

This describes the planned internal migration. Several engineering functions
currently still use dedicated assembly kernels; this documentation does not
claim that they already delegate to the DSL. The older implementation can be
replaced while keeping the engineering entry points useful and familiar.
Engineering function names are not classified as deprecated merely because
their current backend is older. Using the DSL directly is an additional choice,
not a requirement to rename an existing engineering workflow.

## Engineering workflows

- [Linear mechanics, stress/strain and mechanical eigen helpers](linear.md)
- [Heat assembly and solution](heat.md)
- [Engineering load/support declarations](boundary-conditions.md)

## Legacy implementation and formulation mapping

The dedicated [Poisson](poisson.md) and [nonlinear assembly](nonlinear.md)
routines are documented as legacy implementations. Shared field calculus and
[energy-based constitutive utilities](constitutive.md) remain part of modern
DSL workflows even when their definitions live in the same source files.

The table relates familiar engineering/legacy operations to their weak forms.
It is useful for understanding the backend or writing a custom model; it is
not a migration requirement for users of the retained engineering API.

| Operation | DSL formulation / topic |
|:--|:--|
| `poissonMatrix(P)` | `∫(Grad(P) ⋅ Grad(P))` for scalar diffusion |
| `reactionMatrix(P; coefficient=c)` | `∫(P ⋅ c ⋅ P)` |
| `advectionMatrix(P; dir=1, coefficient=b)` | `∫(P ⋅ [b 0.0] ⋅ Grad(P))` for a 2D scalar unknown |
| `gradDivMatrix(U)` | `∫(Div(U) ⋅ Div(U))` |
| `curlCurlMatrix(U)` | `∫(Curl(U) ⋅ Curl(U))` |
| Elastic stiffness | `∫(SymGrad(U) ⋅ D(:Solid,mat) ⋅ SymGrad(U))` for a 3D unknown |
| `loadVector(P, ...)` | `∫(P ⋅ source)` or `∫(P ⋅ traction; Γ="right")` |
| Static displacement or temperature solve | `solveField(K,f; support=bc)` after assembly |
| Mechanical eigen helpers | `solveEigen(A,B; ...)`, with explicit eigenvalue interpretation |
| Dedicated nonlinear assembly | Constitutive response plus residual/tangent assembly and an equilibrium iteration |

Preserve boundary terms, coefficient dimensions, physical groups and sign
conventions when replacing an implementation. Axisymmetric weights/kinematics,
plane-stress reductions, result types and engineering defaults need particular
care. Shared DSL assembly supports direct CSC storage and threading; actual
performance depends on model, quadrature and thread count.

The [older example gallery](../tutorials/legacy-examples.md) is retained.
Existing transient reference material remains available; this update adds no
new transient workflow.
