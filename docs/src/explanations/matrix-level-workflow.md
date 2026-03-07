# Matrix-Level Workflow

A core LowLevelFEM pattern is explicit access to assembled matrices and vectors,
enabling direct customization of algebraic forms.

## Typical sequence

1. Assemble (`stiffnessMatrix`, `massMatrix`, `loadVector`, ...)
2. Apply constraints (`applyBoundaryConditions!`, `freeDoFs`, ...)
3. Solve (`solveDisplacement`, eigen/time integrators)
4. Postprocess (`solveStress`, `show*Results`, `probe`)

## Status

Placeholder page.
