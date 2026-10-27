# Legacy Nonlinear Assembly API

```@meta
CurrentModule = LowLevelFEM
```

Dedicated total-Lagrangian assembly, stress and follower-load routines retained
for compatibility. Modern [constitutive response and tensor conversion](constitutive.md)
and [DSL nonlinear formulations](../tutorials/nonlinear-large-deformation.md)
are documented separately.

## Kinematics and Constitutive Terms

```@docs
nodePositionVector
materialTangentMatrix
initialStressMatrix
IIPiolaKirchhoffStress
```

## Nonlinear Assembly

```@docs
externalTangentFollower
internalForceVector
```

## Visualization

```@docs
showDeformationResults
```

## Legacy derivative helpers

```@docs
∇_old
grad_xy
```
