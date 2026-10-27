# Matrix-Level Workflow

An assembled matrix carries both numerical entries and the field metadata
needed to interpret those entries. This permits inspection and customization
between assembly, constraint handling and solution.

## Single field

The [scalar diffusion example](../manual/assembly-and-solvers.md) explicitly
constructs `K`, `f` and supports. `K.A` is the original sparse matrix and `f.a`
the assembled nodal RHS. The constrained solver reconstructs a nodal field
on the original DOF space, so constrained rows can carry reactions rather than
satisfying the unconstrained RHS equation.

## Two analytically verifiable fields

This model couples scalar fields through `[2 1; 1 3]` at every point of a line.
The exact constant solution is `(p,q)=(1,2)` because `[2 1;1 3]*[1;2]=[4;7]`.
No gradient or essential constraint is needed for this positive-definite mass
system.

```julia
# docs-check: coupled-scalars
using LowLevelFEM, LinearAlgebra

gmsh.initialize()
try
    line_mesh(n=2, order=1)
    mat = Material("body")
    P = Field([mat], type=:ScalarField, dim=1, fieldName=:p)
    Q = Field([mat], type=:ScalarField, dim=1, fieldName=:q)
    A = ∫(P ⋅ 2.0 ⋅ P)
    B = ∫(P ⋅ Q)
    C = ∫(Q ⋅ 3.0 ⋅ Q)
    K = SystemMatrix([A B; B' C])
    F = SystemVector([∫(P ⋅ 4.0), ∫(Q ⋅ 7.0)])
    p, q = solveField(Symmetric(K), F)
    @assert maximum(abs.(p.a .- 1.0)) < 1e-10
    @assert maximum(abs.(q.a .- 2.0)) < 1e-10
    @assert norm(K.A*vcat(p.a,q.a)-F.a) < 1e-10
finally
    gmsh.finalize()
end
```

Reordering fields requires reordering both matrix blocks and RHS fields.
Results follow that order. A system containing only one block returns its
single field directly.

## Eigenproblems reuse the same assembly

For a scalar line field, stiffness and mass can be assembled independently:

```julia
A = ∫(Grad(P) ⋅ Grad(P))
M = ∫(P ⋅ P)
bc = [BoundaryCondition("left", field=P, p=0.0),
      BoundaryCondition("right", field=P, p=0.0)]
modes = solveEigen(A, M; n=1, solver=:direct, support=bc)
λ1 = eigenValue(modes)
ϕ1 = eigenVector(modes)
```

This continuation requires the line model and descriptor to remain active,
so run it before the `finally` above. With two equal linear elements on the
unit interval, the only free mode has the discrete eigenvalue `12`, rather
than the continuum value `π²`. Mesh refinement is a separate accuracy question.

See [Solution API](../reference/solution.md) for normalization, complex output,
reduced/spectral interpolation, optional backends, local coordinates and MPCs.
