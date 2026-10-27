# Weak-Form DSL Design

The DSL separates field descriptors, differential operators, coefficients,
integration and solution. A descriptor such as `P` identifies an unknown and
its interpolation; `Grad(P)` describes how its element DOFs are differentiated.
No solution values are evaluated until a field is supplied after assembly.

## From weak form to matrix

Scalar diffusion has the bilinear form

```math
a(v,u) = \int_\Omega (\nabla v)^T C\,\nabla u\,\mathrm d\Omega.
```

It is represented by `∫(Grad(P) ⋅ C ⋅ Grad(P); Ω="body")`.
The left operator defines the test/row space, the right defines the trial/column
space. For different descriptors, `∫(Q ⋅ Div(U))` is a rectangular coupling
block. `SystemMatrix` combines compatible blocks into a global system.

At each Gauss point the element kernel forms test/trial operator matrices,
evaluates the coefficient, and accumulates the local bilinear contribution.
The standard path scatters contributions directly into CSC sparse storage.
A linear form such as `∫(P ⋅ g)` instead assembles an equivalent nodal RHS;
it is not merely an array of sampled source values.

## Coefficients and compound operators

Coefficient chains can contain numbers, compatible scalar fields and matrices
of scalar entries/fields. Dimensions must match operator output dimensions;
not every field representation is accepted by every integration path.

Compound sums such as `B = Grad(P) + [1.0; 0.0;;] ⋅ P` represent a differential
operator plus a zeroth-order term on a 2D scalar field. `∫(B' ⋅ B)` expands
into ordinary bilinear terms and reuses structural patterns by field pair.

`full(term)` and `reduced(term)` annotate quadrature in compound forms;
`gauss=:auto` honors these rules and uses full quadrature for mixed-rule
products. This does not change interpolation order (`reducedOrder=true`).

## Boundaries and solution

`Ω` selects a domain group and `Γ` a boundary group; with neither, material
groups define the domain. These selectors are mutually exclusive. Essential
conditions enter `solveField`; natural and Robin terms enter the assembled
weak form. The solver preserves matrix/RHS visibility and reconstructs the
original field coordinates after constraint transformations.

See [Weak-form reference](../reference/multifield.md),
[Assembly and solvers](../manual/assembly-and-solvers.md), and
[Matrix-level workflow](matrix-level-workflow.md).
