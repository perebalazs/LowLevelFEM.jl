# LowLevelFEM – API logical structure (draft, function-level)

⚠️ Working document.
This is NOT documentation, but a planning map for later API restructuring.

---

## Core

Fundamental data structures and concepts.

### Problem and materials

- Problem
- Material
- material
- elastic_constants
- elasticSupport
- elasticSupportMatrix
- printElasticConstantsTable

### Fields

- ScalarField
- VectorField
- TensorField
- scalarField
- vectorField
- tensorField
- field
- mergeFields
- projectScalarField
- fieldError
- loadField
- saveField
- isNodal
- isElementwise
- isSaved

### Mesh and geometry

- Geometry
- CoordinateSystem
- generateMesh
- rotateNodes
- projectTo2D
- expandTo3D
- getDimForPhysicalName
- getTagForPhysicalName

### Mapping and integration

- nodesToElements
- elementsToNodes
- integrate
- normalVector
- normalVector1D
- tangentVector

---

## Preprocessing

Preparation of FEM problems.

### Boundary conditions

- BoundaryCondition
- BoundaryConditionFields
- displacementConstraint
- temperatureConstraint
- pressureConstraint
- elasticSupport
- applyBoundaryConditions
- applyBoundaryConditions!
- applyElasticSupport!
- suppressDeformationAtBoundaries
- suppressDeformationAtBoundaries!

### Initial conditions

- initialDisplacement
- initialDisplacement!
- initialVelocity
- initialVelocity!
- initialTemperature
- initialTemperature!

### Loads and sources

- loadVector
- loadVectorOld
- loadVectorNew
- nonFollowerLoadVector
- equivalentNodalForce
- thermalLoadVector
- thermalLoadVectorAXI
- thermalLoadVectorSolid
- heatSource
- heatSourceVector
- heatFlux
- heatFluxVector
- heatConvection
- heatConvectionVector

---

## Operators

Weak-form building blocks, physics-independent.

### Scalar field operators (Poisson / Laplace)

- poissonMatrix
- poissonMatrixVector
- poissonMatrixSymGradF
- stiffnessMatrixPoissonAllInOne
- traceLaplaceMatrix
- tensorLaplaceMatrix
- reactionMatrix
- sourceVector
- solveField

### Vector field operators

- gradMatrix
- gradDivMatrix
- gradDivMatrixF
- curlCurlMatrix
- symmetricGradientMatrix
- tensorDivDivMatrix
- beltramiMichellMatrix

### Advection and flow-related operators

- advectionMatrix
- navierStokesAdvectionMatrix

### Mass and auxiliary operators

- massMatrix
- mapScalarField
- unitTensor
- estimateLengthOfIJV

---

## Mechanics

### Linear mechanics

Small-deformation solid mechanics.

#### Stiffness and mass

- stiffnessMatrix
- stiffnessMatrixSolid
- stiffnessMatrixSolidParallel
- stiffnessMatrixSolidParallelNew
- stiffnessMatrixAXI
- stiffnessMatrixTruss
- massMatrixSolid
- massMatrixSolidParallel
- massMatrixTruss
- dampingMatrix

#### Constitutive relations

- constitutive_matrix
- constitutive_matrix_stress

#### Solvers

- solveDisplacement
- solveStress
- solveStressSlow
- solveStrain
- solveStrainNew
- solveAxialForce

#### Eigenvalue and stability

- solveEigenModes
- solveModalAnalysis
- solveBuckling
- solveBucklingModes
- largestEigenValue
- smallestEigenValue
- largestPeriodTime
- smallestPeriodTime

---

### Nonlinear mechanics

Large deformation and energy-based formulations.

#### Kinematics

- deformationGradient
- nodePositionVector
- grad
- grad_xy
- div
- curl
- rot

#### Internal forces and stresses

- internalForceVector
- internalForceTL
- stress_from_energy
- IIPiolaKirchhoff

#### Tangent operators

- materialTangentMatrix
- tangentMatrixConstitutive
- tangentMatrixInitialStress
- geometricTangentMatrix
- tangent_from_energy
- externalTangentFollower
- externalTangentFollowerTL

#### Solvers and postprocessing

- solveDeformation
- showDeformationResults

---

### Structural dynamics

Mechanical systems with inertia.

- massMatrix
- dynamicResidual
- dynamicTangent
- HHT
- CDM
- initialStressMatrix

---

## Heat transfer

### Steady-state heat conduction

- heatConductionMatrix
- heatCondMatrixSolid
- heatCondMatrixAXI
- heatFluxVector
- solveHeatFlux

### Transient heat conduction

- heatCapacityMatrix
- latentHeatMatrix
- solveTemperature
- q_free
- q_theta
- FDM

---

## Postprocessing

Evaluation and visualization of results.

### Field visualization

- showScalarResults
- showStressResults
- showStrainResults
- showHeatFluxResults
- showModalResults
- showBucklingResults
- showElementResults
- showOnSurface
- plotOnPath
- probe

### Resultants and derived values

- resultant
- resultant2
- flowRate
- fieldsToVolume

---

## Utilities

Internal helpers, debugging, and infrastructure.

- debugMatrix
- debugVector
- checkConsistency
- make_workspace
- make_ws
- make_ws_axi
- inv2x2_fast
- inv3x3
- inv3x3_fast
- openPreProcessor
- openPostProcessor
- showGapThickness
