# Fields and Results

The DSL workflow returns ordinary scalar, vector or tensor fields:

1. Assemble a matrix and RHS with `∫`, then call `solveField`.
2. Derive quantities with field algebra and `∇`.
3. Inspect components, integrate, or create Gmsh views.

For the 3D example in [First Problem](first-problem.md), `u` is nodal, while
`εtensor` and `σ` are element-wise derived fields. `u[1]` is the first
component field; `σ[1,2]` is tensor shear stress. `showDoFResults(u)` and
`showStressResults(σ)` create view tags for visualization and probing.
Keep the Gmsh model active while evaluating or displaying these fields.

Use `isNodal`, `isElementwise`, `nodesToElements` and `elementsToNodes` to
inspect or convert representation. Smoothing discontinuous element results
onto shared nodes can change them; see [Nodal vs Element Fields](../explanations/nodal-vs-element-fields.md).

Eigenresults use the same containers, but columns hold modes and `.t` stores
eigenvalues. Use `eigenValue`, `eigenVector` and `showEigenResults`; their
columns are not a physical time series.

See [Field algebra](../reference/operators.md), [Field calculus](../reference/field-calculus.md),
and [Postprocessing](../manual/postprocessing-and-visualization.md).
