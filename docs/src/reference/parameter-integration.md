# [Parameter Integration](@id parameter-integration)

```@meta
CurrentModule = LowLevelFEM
```

Use `integrate(f, a, b; order=5)` or `∫(f, a, b; order=5)` to integrate
with respect to an independent scalar parameter. For example, a shell's
through-thickness integral eliminates the thickness coordinate while keeping
the spatial coordinates of the shell midsurface:

```math
R(\mathbf x)=\int_a^b f(\mathbf x,z)\,\mathrm dz.
```

This differs from `integrate(P, "body", f)`, which integrates over the mesh
and returns a number, and from `integrate(field)`, which integrates along the
stored time history. The new overload preserves all spatial field values and
time columns in the result. It uses Gauss–Legendre quadrature with
`order ÷ 2 + 1` points; polynomials of degree up to `order` are integrated
exactly up to roundoff. The default `order=5` uses three points. Integrate
separate intervals when the integrand is discontinuous, for example at layer
interfaces.

## Supported results

| `f(z)` returns | Integral returns |
|:--|:--|
| Number | Number |
| `ScalarField`, `VectorField`, `TensorField` | Same field type in elementwise storage |
| Numeric vector or matrix | Numeric vector or matrix |
| Vector or matrix of numbers and `ScalarField`s | Vector or matrix with the same shape; field entries are elementwise |

An array may mix numeric constants and field entries. Numeric entries remain
numbers unless the corresponding integrand entry becomes a field at another
quadrature point. No mesh is required for purely numeric integration. The
input data are not modified. Reversing the endpoints changes the sign;
identical endpoints return a zero with the integrand's shape.

```julia
I = ∫(z -> z^2, -0.5, 0.5)  # 1/12
q = scalarField(P, "body", (x,y,z)->1+x+y)
Q = [q 0; 0 2.0]
A = ∫(z -> Q, -0.5, 0.5)
D_b = ∫(z -> z^2*Q, -0.5, 0.5)
```

The nonconstant entries of `A` and `D_b` are `ScalarField`s, not a new matrix
field type. A constitutive matrix can have any compatible size, including
6×6. Ordinary `TensorField`s still represent 3×3 tensor values.

All field-valued samples of a component must agree on the Problem, field
layout, element support, and stored time grid. Element tags may occur in a
different order. By default, nodal data are copied to the Problem's material
groups; existing elementwise data retain their element support. For a specific
surface or other group, use `onPhysicalGroup`:

```julia
A = ∫(z -> Q, -0.5, 0.5; onPhysicalGroup="shell")
```

Every participating elementwise component must cover that group. A missing
element is an error. This operation never averages an elementwise field to
extend its domain.

## [Through-thickness integration of a layered shell](@id layered-shell-parameter-integration)

Let `z` contain the layer interfaces, with `z[1] = -h/2` and
`z[end] = h/2`. Let `Q[k]` be the plane-stress constitutive matrix of layer
`k`, already expressed in the common shell coordinate system. For a shell
kinematics with in-plane strain `ε(z) = ε₀ + zκ`, the membrane, coupling and
bending matrices are

```math
A=\sum_k\int_{z_k}^{z_{k+1}}Q_k(z)\,\mathrm dz,\qquad
B=\sum_k\int_{z_k}^{z_{k+1}}zQ_k(z)\,\mathrm dz,\qquad
D=\sum_k\int_{z_k}^{z_{k+1}}z^2Q_k(z)\,\mathrm dz.
```

For constant-through-thickness properties within each layer:

```julia
# Q[k] may contain numbers and spatial ScalarFields on the shell.
A = sum(∫(ζ -> Q[k],     z[k], z[k+1]; onPhysicalGroup="shell") for k in eachindex(Q))
B = sum(∫(ζ -> ζ*Q[k],   z[k], z[k+1]; onPhysicalGroup="shell") for k in eachindex(Q))
D_b = sum(∫(ζ -> ζ^2*Q[k], z[k], z[k+1]; onPhysicalGroup="shell") for k in eachindex(Q))
D_mb = [A B; B D_b]

# Existing compound operator with membrane-strain and curvature rows.
K = ∫(B_mb' ⋅ D_mb ⋅ B_mb; Γ="shell")
```

Equivalently, form the entire block matrix in one integral per layer:

```julia
D_mb = sum(
    ∫(ζ -> [Q[k] ζ*Q[k]; ζ*Q[k] ζ^2*Q[k]], z[k], z[k+1];
      onPhysicalGroup="shell") for k in eachindex(Q)
)
```

Layer-wise integration retains membrane–bending coupling for unsymmetric
layups. For one homogeneous layer, `A = h*Q`, `B = 0` and `D_b = h^3*Q/12`.
Use `Q[k](ζ)` instead of a constant matrix if a layer varies continuously
through its thickness. The example uses the linear thickness kinematics and
flat-section measure; more general kinematics or Jacobian factors belong in
the integrand. Transverse shear and drilling terms may be assembled separately
as in the existing shell workflow.

## Using elementwise coefficients in weak forms

The ordinary weak-form assembler accepts nodal and elementwise scalar
coefficients, including entries of a matrix chain. Elementwise data must cover
every element being integrated. If the field is defined only on part of the
mesh, integrate that physical group, provide explicit zero values on other
elements, or explicitly transfer the field with `elementsToElements` or
`elementsToNodes` followed by `nodesToElements`, as appropriate. The assembler
does not silently smooth or extrapolate missing data. General weak-form sums
accept elementwise coefficients too. The separate `ContactGap` coefficient
interface is unchanged.

See [`integrate`](@ref) for the complete function reference and
[Weak-form assembly](multifield.md) for matrix chains.
