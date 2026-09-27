# Contact Mechanics

This tutorial collection demonstrates frictionless slave–master contact in LowLevelFEM using the `ContactGap` weak-form operator.

The examples cover two-dimensional Hertz-type contact with both penalty and Lagrange-multiplier formulations, as well as a three-dimensional penalty-contact problem. Contact detection is based on closest-point projection in the current configuration, while the contact contribution is assembled directly from the weak form at slave-side Gauss points.

The penalty examples show how the contact tangent can be written compactly as `∫(G ⋅ cn ⋅ G)`, with `G = ContactGap(C)`, and solved incrementally using prescribed displacement steps. The Lagrange-multiplier example demonstrates the corresponding mixed formulation and the use of a reduced-order multiplier field for stable contact-pressure approximation with higher-order displacement elements.

## Examples

[2D Hertz contact — penalty method](https://github.com/perebalazs/LowLevelFEM.jl/blob/main/examples/2D-Hertz-penalty-minimal.ipynb)

[2D Hertz contact — Lagrange multiplier method](https://github.com/perebalazs/LowLevelFEM.jl/blob/main/examples/2D-Hertz-lagrange-minimal.ipynb)

[3D contact — penalty method](https://github.com/perebalazs/LowLevelFEM.jl/blob/main/examples/contact-3D-penalty-minimal.ipynb)

## Related

- [Reference: Contact](../reference/contact.md)
- [Reference: Multifield](../reference/multifield.md)
- [Explanations: Weak-Form DSL](../explanations/weak-form-dsl-design.md)
