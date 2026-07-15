# Postprocessing and Visualization

LowLevelFEM provides built-in postprocessing and visualization tools based on Gmsh. Results can be displayed interactively, exported, probed at arbitrary locations, or evaluated along user-defined paths.

Typical workflow:

```julia
u = solveDisplacement(problem, load=[force], support=[bc])
S = solveStress(u)

showDoFResults(u)
showStressResults(S)

openPostProcessor()
```

The most commonly used visualization functions are:

- `showDoFResults` – displacement and other nodal degrees of freedom
- `showStressResults` – stress components and equivalent stress
- `showStrainResults` – strain components
- `showHeatFluxResults` – heat flux fields
- `showElementResults` – element-wise results (e.g. stress, strain, heat flux), which may be discontinuous across element boundaries
- `plotOnPath` – evaluate and plot results along a user-defined path
- `probe` – evaluate field values at arbitrary locations

For complete function documentation, see the API reference below.

## Related API

- [Reference: Postprocessing](../reference/postprocessing.md)
- [Reference: Extra](../reference/extra.md)
