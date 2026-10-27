# Scalar Diffusion, Reaction and Advection with the DSL

For new scalar PDEs use `Field` plus weak-form assembly. The dedicated
[Poisson matrices](../reference/poisson.md) are legacy compatibility routines.

A complete diffusion model with exact solution `u=x` is provided in
[Assembly and solvers](../manual/assembly-and-solvers.md). Starting from its
2D scalar descriptor `P`, additional bilinear terms are explicit:

```julia
K_diffusion = ∫(Grad(P) ⋅ 2.0 ⋅ Grad(P))
K_reaction = ∫(P ⋅ 3.0 ⋅ P)
K_advection = ∫(P ⋅ [1.0 0.0] ⋅ Grad(P))
f = ∫(P ⋅ 1.0)
```

The row `[1.0 0.0]` contracts the trial gradient with a prescribed velocity.
The current `Adv(P)` convenience constructor is incomplete; use this coefficient
chain instead. Match diffusion tensors and coefficient dimensions to the
spatial dimension, and retain your PDE's integration-by-parts sign convention.

Assemble `K = K_diffusion + K_reaction + K_advection` and apply suitable
essential/natural conditions through the same workflow. Advection can make the
matrix nonsymmetric; this expression introduces no stabilization automatically.
For advection-dominated problems, discretization/stabilization must be chosen
and validated explicitly.

See [Weak forms](../reference/multifield.md), [Field calculus](../reference/field-calculus.md),
and [Legacy migration](../reference/legacy.md).
