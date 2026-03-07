# Linear API

Linear mechanics, dynamic solvers, and modal/buckling workflows.

## Assembly

```@docs
stiffnessMatrix
nonLinearStiffnessMatrix
massMatrix
dampingMatrix
elasticSupportMatrix
loadVector
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

## Eigen and Buckling

```@docs
solveEigenModes
solveBucklingModes
solveModalAnalysis
solveBuckling
largestPeriodTime
smallestPeriodTime
largestEigenValue
smallestEigenValue
getEigenVectors
getEigenValues
```

## Time Integration

```@docs
CDM
HHT
CDMaccuracyAnalysis
HHTaccuracyAnalysis
initialDisplacement
initialDisplacement!
initialVelocity
initialVelocity!
```

## Utility

```@docs
printElasticConstantsTable
```
