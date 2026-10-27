# Field Algebra API

```@meta
CurrentModule = LowLevelFEM
```

These operations evaluate existing field values. To construct an element matrix
from unknown fields, use the [weak-form DSL](multifield.md).

## Arithmetic and contractions

| Expression | Meaning |
|:--|:--|
| `s1 + s2`, `s1*s2`, `s1/s2`, `s^p` | Scalar-field arithmetic |
| `s + a`, `v*a`, `T/a` | Arithmetic with scalar numbers |
| `v1 + v2`, `v*s`, `v/s` | Vector addition and scalar scaling |
| `dot(v,w)` or `v ⋅ w` | Vector inner product, returning a scalar field |
| `v × w` | Three-component cross product |
| `v ∘ w` | Dyadic product, returning a tensor field |
| `T*U` | Matrix product at field sample points |
| `T ⋅ U` | Tensor double contraction, returning a scalar field |
| `T*v`, `v*T` | Tensor-vector products |
| `transpose(T)`, `T'` | Tensor transpose |
| `norm(v)` | Pointwise vector norm |
| `det(T)`, `trace(T)`, `inv(T)` | Tensor determinant, trace and inverse |
| `sqrt(T)`, `cbrt(T)`, `log(T)` | Small-matrix functions at field sample points |
| `eigen(T)` | Principal values (`VectorField`) and vectors (`TensorField`) |

Most field algebra evaluates element-wise data, converting nodal inputs when
necessary. Operands must have compatible mesh support, element ordering and
columns. Ordinary array broadcasting is not a substitute for field operations.
Matrix functions retain their usual domain restrictions; for example a real
matrix logarithm is not defined for every tensor.

`eigen(T)` computes local tensor eigensystems. The global FEM eigenproblem is
[`solveEigen`](solution.md).

## Vector and tensor operations

```@docs
dot(::VectorField,::VectorField)
∘(::VectorField,::VectorField)
*(::TensorField,::VectorField)
*(::VectorField,::TensorField)
```

## Tensor matrix functions

```@docs
eigen(::TensorField)
sqrt(::TensorField)
cbrt(::TensorField)
log(::TensorField)
```

Three-component vector contractions require `:v3D` data; expand a
two-component nodal field explicitly before `dot`, dyadic or tensor-vector
products. Tensor `sqrt`, `cbrt` and `log` are matrix functions; scalar-field
versions operate on individual scalar values.

## Scalar and tensor helpers

```@docs
⋅
unitTensor
trace
det(::TensorField)
mapScalarField
```

The `ε` name is an alias of the DSL [SymGrad](multifield.md) constructor. It is
not a stored strain field and does not replace tensor shear conventions.

## Tensor component conversion

[tensorToMatrix and tensorToVector](constitutive.md) expose scalar components
for coefficient construction; conversion to Voigt form does not double shear
entries. See [component conventions](multifield.md#Voigt-convention).

## Related topics

- [Field calculus](field-calculus.md)
- [Fields and storage](fields.md)
- [Eigenvalue transformations](solution.md)
