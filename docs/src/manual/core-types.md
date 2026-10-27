# Core Types

## Model descriptors and values

`Material("body"; ...)` associates properties with a named physical group.
`Field([mat]; type=:ScalarField, dim=2, fieldName=:T)` defines an unknown and
its finite-element interpolation on the current Gmsh model. `Field` and
`Problem` name the same descriptor type; neither contains a solved temperature.

`ScalarField`, `VectorField` and `TensorField` contain actual values and refer
back to a descriptor through `.model`. Lowercase constructors such as
`scalarField` create nodal values, while uppercase constructors normally create
values at element nodes. See [Fields and storage](../reference/fields.md).

A 2D vector unknown has two nodal DOFs per node. Element-wise vector quantities
may use three Cartesian components. Tensors store all nine entries of a 3×3
matrix in column-major order even on a 2D mesh. Do not infer vector/tensor
storage dimensions from the mesh dimension alone.

## Algebraic systems

`SystemMatrix` combines a sparse `.A` matrix with trial/test field metadata.
A rectangular coupling block can connect different fields. `SystemMatrix`
can also assemble a square block system whose field order is stored in
`.problems` and `.offsets`.

`SystemVector` concatenates nodal RHS matrices in that same field order.
All columns are retained as independent right-hand sides. `solveField` returns
a nodal solution field for a single unknown or a tuple for multiple unknowns;
a one-block system returns its sole field directly.

## Columns and eigenfields

For ordinary fields, `.t` describes stored steps. For eigenfields,
`.kind` is `:eigen_real` or `:eigen_imag`, columns are modes and `.t` contains
eigenvalue components. Use `eigenValue` and `eigenVector` to query modes;
do not treat them as a time history.

See [Solution API](../reference/solution.md) and
[Nodal vs Element Fields](../explanations/nodal-vs-element-fields.md).
