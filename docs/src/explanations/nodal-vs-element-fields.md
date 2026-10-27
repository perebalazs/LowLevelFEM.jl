# Nodal vs Element Fields

The same finite-element function can be stored globally at mesh nodes or
locally at each element's nodes. The choice matters when results are
discontinuous or only cover a subset of the model.

| Representation | Storage | Meaning at a shared node |
|:--|:--|:--|
| Nodal | `f.a`, rows are global nodal DOFs | One value per global node/component |
| Element-wise | `f.A[e]`, element tags in `f.numElem` | Each neighboring element may have a different value |

Each matrix column is a stored step/RHS; for eigenfields it is a mode. Inspect
`isNodal(f)` and `isElementwise(f)` instead of reading both arrays as populated.

## Construction and conversion

For a scalar descriptor `P` on a mesh containing `"body"`:

```julia
sn = scalarField(P, "body", (x,y,z) -> x + 2y)
se = nodesToElements(sn)
@assert isNodal(sn)
@assert isElementwise(se)
smooth = elementsToNodes(se)
```

Lowercase construction is nodal; uppercase `ScalarField(P,"body",...)` is
normally element-wise. `nodesToElements` evaluates the global nodal field on
each element without inventing jumps. `elementsToNodes` combines values from
adjacent elements according to its projection/smoothing rule; it can lose jumps.
A discontinuous field is therefore not generally recovered by a round trip.

## Why derivatives are element-wise

A continuous piecewise-linear scalar field has a constant gradient inside each
element, but gradients can differ across an interface. `∇(s)` preserves these
separate values. An element-wise input is differentiated directly, including
its stored support; a field restricted to a surface remains on that surface.
Smoothing before differentiation changes the field being differentiated.

## Components and support

2D vector unknowns have two nodal DOFs per node, while element-wise Cartesian
vector data may use three components. Full tensors store nine components
column-major, `[xx,yx,zx,xy,yy,zy,xz,yz,zz]`. Use component indexing and public
conversion routines rather than assuming a reshape convention.

`elementsToElements` transfers/restricts fields to specified physical groups;
it is separate from changing between global nodal and element storage.
Preserve compatible element tags, ordering and column counts when combining
fields or using them as coefficients.

See [Fields API](../reference/fields.md), [Field calculus](../reference/field-calculus.md),
and [component ordering](../reference/multifield.md).
