# Space–Time FEM for Elastic Wave Propagation

This tutorial demonstrates a space–time finite element formulation for one-dimensional elastic wave propagation. Instead of discretizing space first and advancing the solution with a time-stepping algorithm, time is treated as an additional finite-element coordinate and the complete impact event is obtained from a single coupled finite-element solve.

The example considers the impact of an elastic bar against a rigid wall. The longitudinal elastodynamic equations are written as a first-order system in particle velocity and normalized axial stress, then discretized over the two-dimensional space–time domain using a least-squares weak formulation and standard continuous Lagrange finite elements.

The solution captures the compression wave, its reflection from the traction-free end, the release wave, and the subsequent rebound of the bar. The example also illustrates a coupled two-field formulation, block-matrix assembly from bilinear forms, initial conditions imposed on the space–time boundary, and characteristic wave propagation without a dedicated transient solver or time-stepping loop.

## Example

[space-time-FEM.ipynb](https://github.com/perebalazs/LowLevelFEM.jl/blob/main/examples/space-time-FEM.ipynb)

## Related

- [Transient Elasticity (2D)](transient-elasticity-2d.md)
- [Weak-Form DSL: Navier-Stokes](multifield-weak-form-dsl.md)
- [Reference: Operators](../reference/operators.md)
