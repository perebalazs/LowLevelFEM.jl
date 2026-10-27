# Field Calculus API

```@meta
CurrentModule = LowLevelFEM
```

Differentiate computed field values. These functions return fields; the
capitalized `Grad`, `Div`, `Curl`, and `SymGrad` constructors instead describe
operators inside an [assembly expression](multifield.md).

## Spatial derivatives

```@docs
∇
∂x
∂y
∂z
grad
curl
rot
div(::Union{VectorField,TensorField})
```

The symbolic notation `u ∘ ∇`, `∇ ∘ u`, `∇ ⋅ u`, and `∇ × u` is implemented
for supported field layouts. The order of a dyadic gradient determines whether
the tensor gradient is transposed; see the [calculus guide](../manual/operators-and-fields.md).
Derivatives belong to elements and can jump across shared nodes.

## Derivatives along stored time columns

```@docs
∂t
```

For eigenfields, columns represent modes and `t` stores eigenvalues. Treating
these columns as a physical time history is invalid.

## Integrals

Use [integrate](postprocessing.md) for scalar spatial integration over a group.
`integrate(field)` and `∫(field)` are time-history integrals; `∫(weak_form)`
assembles matrices or loads. These overloads have different meanings.

## Compatibility helpers

The old `∇_old` and `grad_xy` routines remain in the
[legacy nonlinear reference](nonlinear.md). New postprocessing should use
`∇` and the directional wrappers above.
