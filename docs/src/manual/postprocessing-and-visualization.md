# Postprocessing and Visualization

## Derive fields from a DSL solution

For a 2D displacement field, expand the nodal components before taking a full
3×3 tensor gradient:

```julia
u3 = expandTo3D(u)
εtensor = (u3 ∘ ∇ + ∇ ∘ u3) / 2
```

For a 3D isotropic solid (or plane strain after this expansion), the stress is

```julia
σ = mat.E/(1+mat.ν) *
    (εtensor + mat.ν/(1-2mat.ν)*trace(εtensor)*unitTensor(εtensor))
```

This is not the plane-stress constitutive law. Material laws, tensor/engineering
shear and Voigt order must match the formulation. For scalar conduction use
`q = -k*∇(T)` on a compatible scalar temperature field.

## Inspect and integrate

Use component fields (`u[1]`, `σ[1,2]`) for algebra and derived quantities.
`integrate(P,"body",s)` computes a scalar spatial integral; `resultant` sums
assembled nodal vector entries and is not a general surface quadrature routine.
`probe` queries a Gmsh view: first create a view with a `show*` function.
See individual signatures in [Postprocessing API](../reference/postprocessing.md).

Derivative/stress fields are usually element-wise and may be discontinuous.
`elementsToNodes` smooths shared-node values; use it deliberately when a nodal
representation is needed, preserving the original element field when jumps matter.

## Views and interactive tools

`showDoFResults` creates nodal result views, while `showElementResults` handles
element fields. `showStressResults`, `showStrainResults` and `showHeatFluxResults`
provide specialized components and invariants. View creation returns a Gmsh view
tag; `plotOnPath`, `showOnSurface`, `plotOnBeam` and `probe` use those views as
described in their docstrings. `openPostProcessor()` launches the interactive GUI.

A numerical script can omit GUI launch while still constructing and checking
fields and views. Use [saveField/loadField](../reference/fields.md) for field
persistence, and [rotateResult](../reference/boundary-conditions.md) for local
result components. Eigenmode visualization has its own
[showEigenResults](../reference/solution.md) interface.

The engineering `solveStress`/`solveStrain`/`solveHeatFlux` entry points remain
available and are documented under their [engineering API](../reference/legacy.md)
topics. Their implementation can migrate to the same field/DSL engine.
