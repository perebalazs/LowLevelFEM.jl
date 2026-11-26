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

## Store fewer output samples

Use `save_every` to reduce the returned displacement/velocity history without
changing the time increment or skipping integration steps:

```julia
u, v = HHT(K, M, C, f, u0, v0, n, Δt; support=bc, save_every=10)
# The same keyword is available in CDM and FDM.
```

The initial and final states are always included. Use the returned `.t` labels
for plots, since the last saved interval may be shorter. Sampled loads and
supports must still contain one constant column or all `n` integration samples.
See the [storage and time-grid guide](@ref transient-saving).

## Related

- [Nonlinear Cantilever Dynamics: Secant vs Newton-Type HHT](cantilever-hht-secant-vs-newton.md)

- [Legacy Example Gallery](legacy-examples.md)
- [Reference: Linear](../reference/linear.md)
