# Engineering Linear Mechanics API

```@meta
CurrentModule = LowLevelFEM
```

Problem-oriented assembly and solution functions for the engineering workflow.
Names such as `stiffnessMatrix` and `solveDisplacement` are intended to remain
while their current dedicated kernels are migrated to DSL-backed wrappers.
The current implementation is not yet uniformly DSL-based. For custom forms,
use [weak-form assembly](multifield.md) and [solveField](solution.md) directly.
See [API layers and legacy implementations](legacy.md).

## Assembly

```@docs
stiffnessMatrix
nonLinearStiffnessMatrix
massMatrix
dampingMatrix
elasticSupportMatrix
loadVector
sourceVector
applyBoundaryConditions
applyBoundaryConditions!
applyElasticSupport!
```

## Static and Stress/Strain Solvers

```@docs
solveDisplacement
solveStrain
solveStress
solveAxialForce
```

## Engineering eigenvalue helpers

For the generic solver and eigenfields, use [Solution](solution.md). These
helpers provide specialized mechanical interpretations and the existing
`Eigen` result container for engineering workflows.

```@docs
solveEigenModes
solveEigenFields
solveBucklingModes
solveModalAnalysis
solveBuckling
largestPeriodTime
smallestPeriodTime
largestEigenValue
smallestEigenValue
```

## Time Integration

Start with the [Time Integration guide](transient.md) for method comparison,
all keyword defaults, prescribed histories and [choosing Δt](@ref transient-timestep).
The methods accept general scalar/vector/tensor and multifield operators.

### Central differences (CDM)

```@docs
CDM
```

### HHT / Newmark

```@docs
HHT
```

### Modal accuracy analysis

```@docs
CDMaccuracyAnalysis
HHTaccuracyAnalysis
```

### Initial fields

```@docs
initialDisplacement
initialDisplacement!
initialVelocity
initialVelocity!
```

## Utility

```@docs
printElasticConstantsTable
```
