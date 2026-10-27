# Engineering Heat API

```@meta
CurrentModule = LowLevelFEM
```

Problem-oriented heat matrices, loads and solvers form the engineering API.
Their names are intended to remain while dedicated kernels are migrated to
DSL-backed wrappers. Current implementations still use dedicated heat assembly.
For custom diffusion/heat forms, [Grad and weak-form assembly](multifield.md)
and [solveField](solution.md) are available directly. Existing transient
routines remain documented below; see [API layers](legacy.md).

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
