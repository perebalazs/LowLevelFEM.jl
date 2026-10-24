# Linear API

Linear mechanics, dynamic solvers, and modal/buckling workflows.

## Assembly

```@docs
stiffnessMatrix
nonLinearStiffnessMatrix
massMatrix
dampingMatrix
consistentToLumped
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

`solveEigen` is the generic generalized eigenvalue solver. Its default output is
an eigenfield whose columns are mode shapes and whose `t` vector stores
individual eigenvalues. This is distinct from the physical meaning of a
buckling load multiplier or a frequency, which depends on the matrix pencil.

```@docs
solveEigen
showEigenResults
eigenValue
eigenVector
getEigenValues
LowLevelFEM.mapEigenValues
solveEigenProblem
solveEigenModes
solveBucklingModes
solveModalAnalysis
solveBuckling
largestPeriodTime
smallestPeriodTime
largestEigenValue
smallestEigenValue
```

See [Linear Buckling and Prestressed Natural Frequencies](../tutorials/buckling-and-prestressed-modes.md)
for an example of eigenvalue mappings in two mechanical applications.

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
