# Fields and Results

LowLevelFEM exposes scalar, vector, and tensor field objects directly.

## Common Output Flow

1. Solve primary unknowns (`solveDisplacement`, `solveTemperature`, etc.)
2. Derive secondary fields (`solveStress`, `solveStrain`, `solveHeatFlux`)
3. Visualize and query results (`show*Results`, `plotOnPath`, `probe`)

## Example

```julia
u = solveDisplacement(problem, load=[ld], support=[bc])
S = solveStress(u)

ux = showDoFResults(u, :ux)
s  = showStressResults(S, :s)
```

## See Also

- [Manual: Postprocessing and Visualization](../manual/postprocessing-and-visualization.md)
- [Reference: Postprocessing API](../reference/postprocessing.md)
