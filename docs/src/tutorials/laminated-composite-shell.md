# Laminated Composite Shells

This tutorial extends the displacement–rotation shell formulation to a
five-ply orthotropic laminate. Each ply has its own fiber orientation, and
its constitutive stiffness is rotated into the global shell axes before
integration through the thickness. The rectangular midsurface is generated
with `structured_rect_mesh`; no dedicated composite shell element or
multipoint constraints are required.

## Example

[Laminated composite shell notebook](https://github.com/perebalazs/LowLevelFEM.jl/blob/main/examples/laminated-composite-shell.ipynb)

The 200 × 80 mm cantilever uses five 0.5 mm plies with the stacking sequence
[0°, 45°, 90°, −45°, 0°]. The lamina Young's moduli are 135 GPa along the
fibers and 10 GPa across them. The left edge is clamped, and a total tensile
force of 2000 N and bending moment of 600 N·mm are distributed along the right
edge.

## Through-thickness integration and weak form

The membrane, membrane–bending coupling and bending matrices are obtained
by integrating the rotated plane-stress stiffness of each ply. With `Qlay`
containing the rotated ply matrices and `z` the layer interfaces, the code
follows the laminate stiffness moments directly:

```julia
A_lam = sum(∫(ζ -> Qlay[k],     z[k], z[k+1]) for k in eachindex(Qlay))
B_lam = sum(∫(ζ -> ζ*Qlay[k],   z[k], z[k+1]) for k in eachindex(Qlay))
D_lam = sum(∫(ζ -> ζ^2*Qlay[k], z[k], z[k+1]) for k in eachindex(Qlay))
D_mb = [A_lam B_lam; B_lam D_lam]
```

The unsymmetric stacking sequence produces a nonzero coupling matrix
`B_lam`. The combined membrane and bending stiffness is assembled using the
existing shell operators and matrix chains; transverse shear uses a separate
through-thickness stiffness and reduced quadrature.

## Results and ply stress recovery

The notebook reconstructs the membrane and bending resultants, transverse
shear forces and in-plane stresses in each ply's material axes. In particular,
the fiber-direction stress is evaluated separately for every ply, rather than
averaged through the thickness. The example illustrates the deformation and
stress redistribution due to the stacking sequence without applying a
composite failure criterion.

The saved 40 × 16 element solution gives an average loaded-edge extension
of approximately 0.03024 mm and an out-of-plane displacement of −1.04725 mm,
with rotation about both in-plane axes.

## Related

- [Thin-Walled U-Channel Cantilever (Shell)](u-channel-shell.md) – Shell kinematics, selective quadrature and stress resultant recovery.
- [Parameter Integration](../reference/parameter-integration.md) – Layerwise integration of numbers, matrices and spatial fields.
- [Reference: Multifield](../reference/multifield.md) – Compound weak-form operators and matrix chains.
