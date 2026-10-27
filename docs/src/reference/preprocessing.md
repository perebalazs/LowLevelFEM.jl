# Preprocessing API

```@meta
CurrentModule = LowLevelFEM
```

Set up the current Gmsh mesh, material groups, and unknown-field descriptors.
`Field` is an alias of `Problem`; a descriptor defines interpolation and DOFs,
while `ScalarField`/`VectorField`/`TensorField` contain values.

## Materials and unknowns

```@docs
Material
material
Problem
Field
```

## Mesh and geometry

The structured mesh helpers create named physical groups; see their individual
lists before selecting `Ω`, `Γ`, or a boundary condition.

```@docs
line_mesh
structured_rect_mesh
structured_box_mesh
renumberNodes!
openGeometry
```

## Parameters and interactive geometry

```@docs
setParameter
setParameters
getParameter
openPreProcessor
```

## Next steps

- [Field construction and conversion](fields.md)
- [Boundary conditions, DOFs and coordinates](boundary-conditions.md)
- [Weak-form operators and assembly](multifield.md)
- [Solution and spectral/reduced interpolation](solution.md)
