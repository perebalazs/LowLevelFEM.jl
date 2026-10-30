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
```

## Time Integration

The [Time Integration guide](transient.md) explains FDM inputs, theta choices,
all keyword defaults and [choosing Δt](@ref transient-timestep). The same method
also accepts general first-order scalar/vector/tensor or multifield systems.
Here C denotes heat capacity (A1), while K denotes conduction (A0).

```@docs
FDM
```
