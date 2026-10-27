# Fields and Storage API

```@meta
CurrentModule = LowLevelFEM
```

Lowercase constructors create nodal values; uppercase constructors normally
create element-wise values. Both are finite-element fields associated with a
`Problem`. Inspect the representation rather than inferring it from a variable
name. See [Nodal vs Element Fields](../explanations/nodal-vs-element-fields.md).

## Containers and construction

```@docs
ScalarField
VectorField
TensorField
field
scalarField
vectorField
tensorField
```

## Access, conversion and restriction

`f[i]` selects a vector component; `T[i,j]` selects a tensor component. These
are scalar fields. For column/mode selection use the documented constructor or
[eigenVector](solution.md), not the component indexing syntax.

```@docs
DoFs
isNodal
isElementwise
elementsToNodes
nodesToElements
elementsToElements
mergeFields
projectTo2D
expandTo3D
```

## Persistence

```@docs
saveField
loadField
isSaved
```

## Related topics

- [Field algebra](operators.md): arithmetic, contractions, invariants and tensor functions.
- [Field calculus](field-calculus.md): derivatives of computed fields.
- [Solution](solution.md): `SystemMatrix`, `SystemVector`, `solveField` and `solveEigen`.
