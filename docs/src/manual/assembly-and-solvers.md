# Assembly and Solvers

LowLevelFEM separates finite element assembly from the solution process. Global matrices can be assembled explicitly, inspected or modified, then passed to direct or iterative solvers.

A typical workflow is:

```julia
K = stiffnessMatrix(problem)
f = loadVector(problem, [force])

u = solveDisplacement(K, f, support=[bc])
```

Alternatively, the complete workflow can be executed with a single high-level function:

```julia
u = solveDisplacement(problem,
    load=[force],
    support=[bc])
```

Available analysis types include:

- static linear analysis
- transient dynamics
- modal analysis
- linear buckling
- nonlinear solid mechanics
- coupled multiphysics problems

Direct and iterative linear solvers are supported. Iterative methods can be combined with arbitrary preconditioners provided by Julia's linear algebra ecosystem.

For implementation details and advanced solver options, see the API reference below.

## Related API

- [Reference: Linear](../reference/linear.md)
- [Reference: Nonlinear](../reference/nonlinear.md)
- [Reference: Multifield](../reference/multifield.md)
