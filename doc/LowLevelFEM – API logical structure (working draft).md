# LowLevelFEM – API logical structure (working draft)

## Core

Fundamental data structures and concepts used everywhere.

- Problem
- Material
- BoundaryCondition (base type)
- InitialCondition
- ScalarField
- TensorField
- TimeStepper / TimeGrid
- Mesh / entity access helpers
- Nodes ↔ Elements mapping
  - nodesToElements
  - elementsToNodes

---

## Preprocessing

Preparation of the FEM problem before assembly and solution.

### Mesh handling

- loadMesh
- setCurrentMesh
- entity queries (physical groups, entities)

### Boundary conditions

- BoundaryCondition
- DirichletBC
- NeumannBC
- RobinBC
- FollowerLoadBC
- applyBoundaryConditions

### Initial conditions

- InitialCondition
- setInitialField

### Loads and sources

- BodyForce
- SurfaceLoad
- VolumeSource

---

## Operators

Weak-form building blocks, independent of specific physics.

### Scalar field operators

(typically Poisson-type problems)

- poissonMatrix
- poissonMatrixVector
- laplaceMatrix
- sourceVector
- diffusionOperator

### Vector field operators

- gradDivMatrix
- curlCurlMatrix
- divergenceOperator
- gradientOperator

### Mass and time-related operators

- massMatrix
- dampingMatrix

---

## Mechanics

Solid mechanics formulations.

### Linear mechanics

Small deformation theory.

- stiffnessMatrix
- internalForceVector (linear)
- planeStressMatrix
- planeStrainMatrix
- axisymmetricMatrix
- linearElasticMaterial

### Nonlinear mechanics

Large deformation, energy-based formulations.

- deformationGradient
- greenLagrangeStrain
- firstPiolaKirchhoff
- secondPiolaKirchhoff
- internalForceVector (nonlinear)
- materialTangentMatrix
- geometricTangentMatrix
- followerLoadVector
- followerLoadTangent

### Structural dynamics

Mechanical systems with inertia.

- massMatrix (mechanical)
- dynamicResidual
- dynamicTangent
- NewmarkIntegrator
- HHTIntegrator

---

## Heat transfer

Thermal problems.

### Steady-state heat conduction

- heatStiffnessMatrix
- heatFluxVector
- steadyHeatSolver

### Transient heat conduction

- heatCapacityMatrix
- transientHeatResidual
- transientHeatTangent
- timeStepHeatSolver

---

## Dynamics

General time integration framework (physics-independent).

- TimeIntegrator
- Newmark
- HHT
- generalizedAlpha
- timeResidual
- timeTangent
- advanceTimeStep

---

## Postprocessing

Evaluation and interpretation of results.

### Field evaluation

- evaluateFieldAtNodes
- evaluateFieldAtGaussPoints
- interpolateField

### Derived quantities

- strainField
- stressField
- energyDensity
- heatFlux

### Export and visualization

- exportVTK
- exportCSV
- exportField

---

## Utilities

Internal helpers and debugging tools.

- debugMatrix
- debugVector
- checkConsistency
- estimateLengthOfIJV
- performance helpers
- logging utilities
