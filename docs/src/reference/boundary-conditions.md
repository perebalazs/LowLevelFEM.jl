# Boundary Conditions, MPC and Coordinates

```@meta
CurrentModule = LowLevelFEM
```

For modern DSL models, assemble natural loads with `∫`; pass essential
conditions to `solveField`/`solveEigen`. In a mixed model, bind each condition
to its unknown with `field=P`.

## Conditions and DOF selection

```@docs
BoundaryCondition
BoundaryConditionFields
constrainedDoFs
freeDoFs
allDoFs
ndofs
```

## Multi-point constraints

Detailed construction, rigid links, and periodic examples are in
[Multifield: Multi-point constraints](multifield.md#Multi-point-constraints).

```@docs
MultiPointConstraint
MPC
```

`MPC` is an alias of `MultiPointConstraint`. The rigid-rotation and collapse
helpers are documented with [multifield assembly](multifield.md).

## Local coordinates

```@docs
CoordinateSystem
Transformation
rotateNodes
rotateResult
```

With `coordSys`, ordinary RHS components and prescribed values are interpreted
in the local solver basis; [`Global`](solution.md) marks Cartesian loads.
Solutions are reconstructed in global coordinates. See the
[coordinate guide](../manual/coordinate-systems.md) for a complete example.

## Engineering load and support declarations

These constructors store data for the retained engineering API. Current
engineering assembly routines consume those declarations.
They do not assemble a weak-form load or Robin matrix automatically.

```@docs
LoadCondition
displacementConstraint
load
elasticSupport
temperatureConstraint
heatFlux
heatSource
heatConvection
BoundaryCondition_to_LoadCondition
```
