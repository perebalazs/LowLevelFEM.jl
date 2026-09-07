# Timoshenko Beam Structures

This tutorial demonstrates how planar Timoshenko beam structures can be formulated directly in LowLevelFEM using the weak-form DSL and a coupled displacement–rotation field description.

The examples cover both a general plane beam structure and a Gerber beam with internal hinges. The beam kinematics are defined through a two-component displacement field and an independent cross-section rotation field, while axial deformation, transverse shear, and bending are assembled directly from the weak form. The Gerber beam example additionally demonstrates the use of multi-point constraints (MPCs) to enforce translational continuity at internal hinges while allowing independent rotations.

The tutorials also show how concentrated forces, distributed loads, and concentrated moments can be applied, and how axial force, shear force, and bending moment diagrams can be recovered from the numerical solution.

## Examples

[plane-beam-structure.ipynb](https://github.com/perebalazs/LowLevelFEM.jl/blob/main/examples/plane-beam-structure.ipynb)

[gerber.ipynb](https://github.com/perebalazs/LowLevelFEM.jl/blob/main/examples/Gerber-tutorial.ipynb)

## Related

- [Reference: Multifield](../reference/multifield.md)
- [Explanations: Weak-Form DSL](../explanations/weak-form-dsl-design.md)
