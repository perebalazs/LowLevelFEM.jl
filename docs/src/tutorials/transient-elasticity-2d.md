# Transient Linear Elasticity (2D)

This tutorial demonstrates the transient analysis of a two-dimensional linear elastic solid. In contrast to static analysis, the inertia of the structure is taken into account, allowing the simulation of time-dependent phenomena such as vibrations, wave propagation, and dynamic loading.

The example introduces the complete workflow, including mesh generation, material definition, boundary conditions, time-dependent loading, time integration, and post-processing of displacement and stress fields. 

## Example

[lin_elast_2D_transient.ipynb](https://github.com/perebalazs/LowLevelFEM.jl/blob/main/examples/lin_elast_2D_transient.ipynb)

## Choosing the time increment

CDM must respect the largest frequency of the constrained/reduced system;
lagged damping can tighten that limit. Standard HHT is unconditionally stable
for the usual linear positive-mass system, but both methods still need enough
time samples to resolve the desired response and prescribed motion.

Use the [time-step guide](@ref transient-timestep) for the formulas and a
`largestEigenValue` example. Then halve Δt, rebuild sampled histories and n for
the same final time, and compare displacements, velocities and stress peaks.
The [method and option tables](../reference/transient.md) cover damping,
iterative solves, MPCs and nonlinear operator updates.

## Related

- [Legacy Example Gallery](legacy-examples.md)
- [Reference: Linear](../reference/linear.md)
