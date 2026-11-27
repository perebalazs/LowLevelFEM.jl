# Thin-Walled U-Channel Cantilever (Shell)

This tutorial develops a three-dimensional shell model of a thin-walled
U-channel cantilever directly from its weak form using the LowLevelFEM DSL.
The three connected planar midsurfaces represent the web and flanges; the
unknown displacement and rotation vectors are expressed in global Cartesian
coordinates, with a local orthonormal shell frame on each facet.

## Example

- [U-beam_as_shell.ipynb](https://github.com/perebalazs/LowLevelFEM.jl/blob/main/examples/U-beam_as_shell.ipynb) – Complete tutorial notebook, including the derivation and postprocessing.
- [U-beam_as_shell.geo](https://github.com/perebalazs/LowLevelFEM.jl/blob/main/examples/U-beam_as_shell.geo) – Parametric Gmsh geometry and mesh.

Keep the geometry file in the notebook's working directory. The default
section is 50 mm wide, 100 mm high and 5 mm thick, with a cantilever length of
500 mm. The mesh density can be adjusted in the Gmsh preprocessor.

## Shell formulation

The notebook derives the small-displacement Reissner–Mindlin-type membrane,
bending and transverse-shear operators from the displacement and rotation
fields. The constitutive matrices are obtained by integrating the
homogeneous plane-stress material through the symmetric thickness.
The membrane and bending terms use full quadrature, while the transverse
shear term uses reduced quadrature to alleviate shear locking:

```julia
K = ∫(B_mb' ⋅ D_mb ⋅ B_mb; Γ=surf) +
    ∫(B_s'  ⋅ D3   ⋅ B_s;  Γ=surf, gauss=:reduced)
```

A separate small drilling-rotation stabilization is assembled as `Kd`.
The operators and coefficients are defined in the notebook, without a
predefined U-channel element.

## Rigid end cross-section and results

The clamped root restrains translations and rotations. A remote reference
point at the centroid of the free-end midsurface contour carries a 1000 N
transverse force. Multipoint constraints and `rigidRotationMap` couple the
free-end nodes to the reference point, preserving the rigid-body contribution
of the cross-section rotation. The stiffness and loads are transformed
consistently using `rigidRotationTransformation`.

The notebook recovers the physical displacement, membrane and bending
resultants, transverse shear forces, and in-plane von Mises stresses at the
top, middle and bottom of the shell thickness. The saved solution illustrates
combined bending and torsion: the centroid of an open channel does not
coincide with its shear center, and the rigid end cross-section restrains
warping. The reference-point displacement is approximately 0.2672 mm in the
load direction for the saved mesh; the numerical values depend on mesh
settings and the specified boundary conditions.

## Related

- [Reference: Multifield](../reference/multifield.md) – Compound operators, reduced integration and rigid rotation coupling.
- [Parameter Integration](../reference/parameter-integration.md) – Through-thickness integration for spatially varying or layered sections.
- [Timoshenko Beam Structures](timoshenko-beam-structures.md)
- [Explanations: Weak-Form DSL](../explanations/weak-form-dsl-design.md)
