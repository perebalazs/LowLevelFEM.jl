# Field Algebra, Calculus and DSL Operators

LowLevelFEM supports two related operations: evaluating quantities from existing
fields, and building weak forms for unknown fields. The argument tells you
which operation is taking place.

## Field algebra

Once `u` is a displacement field and `σ` a stress tensor field:

```julia
u3 = expandTo3D(u)  # only when u uses the two-component :v2D layout
εtensor = (u3 ∘ ∇ + ∇ ∘ u3) / 2
energy_density = (σ ⋅ εtensor) / 2
```

`σ ⋅ εtensor` is a double contraction, while `σ*εtensor` is a tensor matrix
product. For vectors, `v ⋅ w` is an inner product and `v ∘ w` is a dyadic
product. See the [field algebra operation table](../reference/operators.md).

## Field calculus

For computed values, `∇(s)` differentiates a scalar field;
`∇(v; nabla=:div)` and `∇(v; nabla=:curl)` evaluate divergence and curl.
Spatial derivatives return element-wise data, preserving jumps between elements.
Supported layouts are listed in [Field calculus](../reference/field-calculus.md).

For a supported vector field, `u ∘ ∇` returns its gradient tensor and
`∇ ∘ u` its transpose. Therefore their symmetric average gives tensor strain.
`∂x(s)`, `∂y(s)` and `∂z(s)` extract Cartesian gradient components. These
are derivatives of the interpolated FE field, not of an unevaluated analytic
function used to construct it.

## Weak-form operators

For an unknown descriptor `U`, `Grad(U)`, `Div(U)`, `Curl(U)`, `SymGrad(U)`
and the surface/axial operators construct symbolic element operators.
The integration call assembles them:

```julia
K = ∫(SymGrad(U) ⋅ D(:PlaneStrain, mat) ⋅ SymGrad(U))
f = ∫(U ⋅ [1.0, 0.0]; Γ="right")
u = solveField(K, f; support=[BoundaryCondition("left", ux=0, uy=0)])
```

Here `U` must describe a 2D vector unknown on an existing mesh.
`SymGrad(U)` uses engineering shear strains; `εtensor` above stores tensor
shears. In 3D the DSL Voigt order is `[xx,yy,zz,xy,yz,zx]`.
Constitutive coefficients must match the chosen convention.

See [weak-form operators and component ordering](../reference/multifield.md),
[constitutive conventions](../reference/constitutive.md), and the
[DSL design explanation](../explanations/weak-form-dsl-design.md).
