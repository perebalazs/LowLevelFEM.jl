# Postprocessing API

```@meta
CurrentModule = LowLevelFEM
```

Visualization, probing, path plots, and result extraction helpers.

## Derived Quantities

```@docs
fieldError
resultant
integrate
normalVector
tangentVector
```

## Result Views

```@docs
showDoFResults
showModalResults
showBucklingResults
showStrainResults
showStressResults
showElementResults
showHeatFluxResults
showOnSurface
plotOnPath
plotOnBeam
```

## Probing and UI

```@docs
probe
openPostProcessor
```

Field persistence is documented under [Fields and storage](fields.md), result
coordinate transformations under [Boundary conditions and coordinates](boundary-conditions.md),
and field derivatives under [Field calculus](field-calculus.md).
