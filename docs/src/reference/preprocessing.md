# Preprocessing API

Preprocessing includes model setup, material/problem definitions,
boundary/load declarations, and mesh-related helpers.

## Mesh and Geometry Setup

```@docs
line_mesh
structured_rect_mesh
structured_box_mesh
projectTo2D
expandTo3D
rotateNodes
```

## Model Setup and Parameters

```@docs
setParameter
setParameters
getParameter
openPreProcessor
openGeometry
```

## Boundary Conditions and Loads

```@docs
displacementConstraint
temperatureConstraint
elasticSupport
load
heatFlux
heatSource
heatConvection
BoundaryConditionFields
BoundaryCondition_to_LoadCondition
```

## DOF Selection Helpers

```@docs
constrainedDoFs
freeDoFs
allDoFs
```
