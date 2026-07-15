# Operators and Fields

LowLevelFEM treats finite element fields as first-class objects. Scalar, vector and tensor fields can be combined with differential operators to construct new quantities directly in Julia code.

Typical examples include:

```julia
ε = (u ∘ ∇ + ∇ ∘ u) / 2
σ = D * ε
q = -k * ∇(T)
```

These expressions closely follow their mathematical notation while remaining fully executable Julia code.

The same operators can also be combined to assemble finite element formulations directly from their weak form. This provides a concise and flexible way to implement custom PDEs and multiphysics problems.

Available functionality includes:

- scalar, vector and tensor finite element fields
- gradient, divergence and curl operators
- field algebra and tensor operations
- field transformations and projections
- weak-form assembly using compound operators

For the complete API and additional examples, see the references below.

## Related API

- [Reference: Operators](../reference/operators.md)
- [Reference: Multifield](../reference/multifield.md)
- [Reference: Fields](../reference/fields.md)
