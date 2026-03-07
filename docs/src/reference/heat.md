# Heat API

Steady/transient heat conduction operators, vectors, and solvers.

## Heat Assembly

```@docs
heatConductionMatrix
heatCapacityMatrix
latentHeatMatrix
heatConvectionMatrix
heatConvectionVector
heatFluxVector
heatSourceVector
thermalLoadVector
applyHeatConvection!
```

## Heat Solvers and Initial Conditions

```@docs
solveTemperature
solveHeatFlux
initialTemperature
initialTemperature!
FDM
```
