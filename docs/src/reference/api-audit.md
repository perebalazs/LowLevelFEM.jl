# Public API Audit

```@meta
CurrentModule = LowLevelFEM
```

This inventory covers the loaded package exports (including aliases, types,
macros and contact constants), not unused historical source files. After this
documentation update, all **224 defined exported bindings** have discoverable
Julia documentation, including inherited alias/macro documentation. Four
additional exported names have no implementation and are listed below.
Docstring presence alone does not prove every behavior is correct; signatures,
keyword variants, examples and component conventions were reviewed separately.

## Known executable-API inconsistencies

These findings are recorded without changing executable `.jl` source:

| Interface | Current behavior / recommended path |
|:--|:--|
| `getEigenVectors` | Exported but undefined; use `eigenVector` for eigenfields, or `output=:raw` from `solveEigen`. |
| `nodalForce!`, `nodalAcceleration!` | Exported but undefined; no callable implementation is supplied by those names. |
| `probe_field` | Exported but undefined; use the documented Gmsh-view `probe` interface where appropriate. |
| `Adv(P)` | Calls the internal constructor without its velocity argument; use an explicit coefficient row times `Grad(P)` as in the [advection guide](../tutorials/poisson-and-custom-operators.md). |
| Energy-based material tangent | The Tensors extension and DSL Voigt orders differ. See the documented two-axis permutation under [constitutive models](constitutive.md). |

`gmsh` is re-exported from the bundled Gmsh API; this project does not duplicate
its complete external reference. `Field`, `MPC`, `ε`, `sourceVector` and
`∫N_c_dΩ` are aliases; their primary documentation describes their target.

## Repeatable docstring check

From the repository root, run:

```sh
python3 docs/check_api_docs.py
python3 docs/check_api_docs.py --strict
```

The default check fails on undocumented defined exports or newly undefined
exports, and reports the four known undefined exports as warnings. `--strict`
also fails for those known undefined names. The checker loads the actual Julia
module and uses `Base.Docs.hasdoc`, so aliases and imported macro documentation
are handled by Julia. It does not edit the project, install packages or run
precompile workloads. It creates a temporary script outside the repository.

The `docs/README.md` file in the source tree describes an isolated
local build and executable-example checks without deployment or dependency
changes. Existing Julia source files are restricted to docstring changes.

## Alphabetical export-to-reference map

| Public name | Primary documentation |
|:--|:--|
| `@disp` | [Reference](extra.md) |
| `@showdef` | [Reference](extra.md) |
| `@showfields` | [Reference](extra.md) |
| `@showmem` | [Reference](extra.md) |
| `@showmethods` | [Reference](extra.md) |
| `@showsize` | [Reference](extra.md) |
| `@showstruct` | [Reference](extra.md) |
| `@showtype` | [Reference](extra.md) |
| `Adv` | [Reference](multifield.md) |
| `AxialGrad` | [Reference](multifield.md) |
| `BoundaryCondition` | [Reference](boundary-conditions.md) |
| `BoundaryConditionFields` | [Reference](boundary-conditions.md) |
| `BoundaryCondition_to_LoadCondition` | [Reference](boundary-conditions.md) |
| `CDM` | [Reference](linear.md) |
| `CDMaccuracyAnalysis` | [Reference](linear.md) |
| `CONTACT_OPEN` | [Reference](contact.md) |
| `CONTACT_SLIP` | [Reference](contact.md) |
| `CONTACT_STICK` | [Reference](contact.md) |
| `ConstitutiveMatrix` | [Reference](constitutive.md) |
| `Contact` | [Reference](contact.md) |
| `ContactGap` | [Reference](contact.md) |
| `ContactSet` | [Reference](contact.md) |
| `ContactStiffness` | [Reference](contact.md) |
| `CoordinateSystem` | [Reference](boundary-conditions.md) |
| `Curl` | [Reference](multifield.md) |
| `D` | [Reference](constitutive.md) |
| `Div` | [Reference](multifield.md) |
| `DoFs` | [Reference](fields.md) |
| `FDM` | [Reference](heat.md) |
| `Field` | [Reference](preprocessing.md) |
| `Global` | [Reference](solution.md) |
| `Grad` | [Reference](multifield.md) |
| `HHT` | [Reference](linear.md) |
| `HHTaccuracyAnalysis` | [Reference](linear.md) |
| `IIPiolaKirchhoff` | [Reference](constitutive.md) |
| `IIPiolaKirchhoffStress` | [Reference](nonlinear.md) |
| `Id` | [Reference](multifield.md) |
| `LoadCondition` | [Reference](boundary-conditions.md) |
| `MPC` | [Reference](boundary-conditions.md) |
| `Material` | [Reference](preprocessing.md) |
| `MultiPointConstraint` | [Reference](boundary-conditions.md) |
| `Problem` | [Reference](preprocessing.md) |
| `ScalarField` | [Reference](fields.md) |
| `SurfaceDiv` | [Reference](multifield.md) |
| `SurfaceGrad` | [Reference](multifield.md) |
| `SurfaceSymGrad` | [Reference](multifield.md) |
| `SymGrad` | [Reference](multifield.md) |
| `SystemMatrix` | [Reference](solution.md) |
| `SystemVector` | [Reference](solution.md) |
| `TangentialGrad` | [Reference](multifield.md) |
| `TensorDiv` | [Reference](multifield.md) |
| `TensorField` | [Reference](fields.md) |
| `Transformation` | [Reference](boundary-conditions.md) |
| `VectorField` | [Reference](fields.md) |
| `advectionMatrix` | [Reference](poisson.md) |
| `allDoFs` | [Reference](boundary-conditions.md) |
| `applyBoundaryConditions` | [Reference](linear.md) |
| `applyBoundaryConditions!` | [Reference](linear.md) |
| `applyElasticSupport!` | [Reference](linear.md) |
| `applyHeatConvection!` | [Reference](heat.md) |
| `beltramiMichellMatrix` | [Reference](poisson.md) |
| `build_csc_pattern` | [Reference](multifield.md) |
| `collapseMPC` | [Reference](multifield.md) |
| `consistentToLumped` | [Reference](multifield.md) |
| `constrainedDoFs` | [Reference](boundary-conditions.md) |
| `contact` | [Reference](contact.md) |
| `curl` | [Reference](field-calculus.md) |
| `curlCurlMatrix` | [Reference](poisson.md) |
| `dampingMatrix` | [Reference](linear.md) |
| `det` | [Reference](operators.md) |
| `displacementConstraint` | [Reference](boundary-conditions.md) |
| `div` | [Reference](field-calculus.md) |
| `eigenValue` | [Reference](solution.md) |
| `eigenVector` | [Reference](solution.md) |
| `elasticSupport` | [Reference](boundary-conditions.md) |
| `elasticSupportMatrix` | [Reference](linear.md) |
| `elementsToElements` | [Reference](fields.md) |
| `elementsToNodes` | [Reference](fields.md) |
| `energyMaterial` | [Reference](constitutive.md) |
| `expandTo3D` | [Reference](fields.md) |
| `externalTangentFollower` | [Reference](nonlinear.md) |
| `field` | [Reference](fields.md) |
| `fieldError` | [Reference](postprocessing.md) |
| `fieldsToVolume` | [Reference](extra.md) |
| `freeDoFs` | [Reference](boundary-conditions.md) |
| `full` | [Reference](multifield.md) |
| `getEigenValues` | [Reference](solution.md) |
| `getEigenVectors` | Exported name; no implementation |
| `getParameter` | [Reference](preprocessing.md) |
| `gmsh` | Bundled Gmsh API (external module) |
| `grad` | [Reference](field-calculus.md) |
| `gradDivMatrix` | [Reference](poisson.md) |
| `gradMatrix` | [Reference](poisson.md) |
| `grad_xy` | [Reference](nonlinear.md) |
| `heatCapacityMatrix` | [Reference](heat.md) |
| `heatConductionMatrix` | [Reference](heat.md) |
| `heatConvection` | [Reference](boundary-conditions.md) |
| `heatConvectionMatrix` | [Reference](heat.md) |
| `heatConvectionVector` | [Reference](heat.md) |
| `heatFlux` | [Reference](boundary-conditions.md) |
| `heatFluxVector` | [Reference](heat.md) |
| `heatSource` | [Reference](boundary-conditions.md) |
| `heatSourceVector` | [Reference](heat.md) |
| `initialDisplacement` | [Reference](linear.md) |
| `initialDisplacement!` | [Reference](linear.md) |
| `initialStressMatrix` | [Reference](nonlinear.md) |
| `initialTemperature` | [Reference](heat.md) |
| `initialTemperature!` | [Reference](heat.md) |
| `initialVelocity` | [Reference](linear.md) |
| `initialVelocity!` | [Reference](linear.md) |
| `integrate` | [Reference](postprocessing.md) |
| `internalForceVector` | [Reference](nonlinear.md) |
| `isElementwise` | [Reference](fields.md) |
| `isNodal` | [Reference](fields.md) |
| `isSaved` | [Reference](fields.md) |
| `largestEigenValue` | [Reference](linear.md) |
| `largestPeriodTime` | [Reference](linear.md) |
| `latentHeatMatrix` | [Reference](heat.md) |
| `line_mesh` | [Reference](preprocessing.md) |
| `load` | [Reference](boundary-conditions.md) |
| `loadField` | [Reference](fields.md) |
| `loadTensor` | [Reference](poisson.md) |
| `loadVector` | [Reference](linear.md) |
| `mapEigenValues` | [Reference](solution.md) |
| `mapScalarField` | [Reference](operators.md) |
| `massMatrix` | [Reference](linear.md) |
| `material` | [Reference](preprocessing.md) |
| `materialTangentMatrix` | [Reference](nonlinear.md) |
| `mergeFields` | [Reference](fields.md) |
| `navierStokesAdvectionMatrix` | [Reference](poisson.md) |
| `ndofs` | [Reference](boundary-conditions.md) |
| `nodalAcceleration!` | Exported name; no implementation |
| `nodalForce!` | Exported name; no implementation |
| `nodePositionVector` | [Reference](nonlinear.md) |
| `nodesToElements` | [Reference](fields.md) |
| `nonLinearStiffnessMatrix` | [Reference](linear.md) |
| `normalVector` | [Reference](postprocessing.md) |
| `openGeometry` | [Reference](preprocessing.md) |
| `openPostProcessor` | [Reference](postprocessing.md) |
| `openPreProcessor` | [Reference](preprocessing.md) |
| `plotOnBeam` | [Reference](postprocessing.md) |
| `plotOnPath` | [Reference](postprocessing.md) |
| `poissonMatrix` | [Reference](poisson.md) |
| `printElasticConstantsTable` | [Reference](linear.md) |
| `probe` | [Reference](postprocessing.md) |
| `probe_field` | Exported name; no implementation |
| `projectScalarField` | [Reference](extra.md) |
| `projectTo2D` | [Reference](fields.md) |
| `reactionMatrix` | [Reference](poisson.md) |
| `reduced` | [Reference](multifield.md) |
| `reduced_bc_data` | [Reference](solution.md) |
| `reductionMatrices` | [Reference](solution.md) |
| `renumberNodes!` | [Reference](preprocessing.md) |
| `resultant` | [Reference](postprocessing.md) |
| `rigidRotationMap` | [Reference](multifield.md) |
| `rot` | [Reference](field-calculus.md) |
| `rotateNodes` | [Reference](boundary-conditions.md) |
| `rotateResult` | [Reference](boundary-conditions.md) |
| `saveField` | [Reference](fields.md) |
| `scalarField` | [Reference](fields.md) |
| `setParameter` | [Reference](preprocessing.md) |
| `setParameters` | [Reference](preprocessing.md) |
| `showBucklingResults` | [Reference](postprocessing.md) |
| `showDeformationResults` | [Reference](nonlinear.md) |
| `showDoFResults` | [Reference](postprocessing.md) |
| `showEigenResults` | [Reference](solution.md) |
| `showElementResults` | [Reference](postprocessing.md) |
| `showHeatFluxResults` | [Reference](postprocessing.md) |
| `showModalResults` | [Reference](postprocessing.md) |
| `showOnSurface` | [Reference](postprocessing.md) |
| `showStrainResults` | [Reference](postprocessing.md) |
| `showStressResults` | [Reference](postprocessing.md) |
| `smallestEigenValue` | [Reference](linear.md) |
| `smallestPeriodTime` | [Reference](linear.md) |
| `solveAxialForce` | [Reference](linear.md) |
| `solveBuckling` | [Reference](linear.md) |
| `solveBucklingModes` | [Reference](linear.md) |
| `solveDisplacement` | [Reference](linear.md) |
| `solveEigen` | [Reference](solution.md) |
| `solveEigenFields` | [Reference](linear.md) |
| `solveEigenModes` | [Reference](linear.md) |
| `solveEigenProblem` | [Reference](solution.md) |
| `solveField` | [Reference](solution.md) |
| `solveHeatFlux` | [Reference](heat.md) |
| `solveModalAnalysis` | [Reference](linear.md) |
| `solveStrain` | [Reference](linear.md) |
| `solveStress` | [Reference](linear.md) |
| `solveTemperature` | [Reference](heat.md) |
| `sourceVector` | [Reference](linear.md) |
| `spectralTransformation` | [Reference](solution.md) |
| `stiffnessMatrix` | [Reference](linear.md) |
| `structured_box_mesh` | [Reference](preprocessing.md) |
| `structured_rect_mesh` | [Reference](preprocessing.md) |
| `subSystemMatrix` | [Reference](solution.md) |
| `symmetricGradientMatrix` | [Reference](poisson.md) |
| `tangentMatrix` | [Reference](constitutive.md) |
| `tangentVector` | [Reference](postprocessing.md) |
| `temperatureConstraint` | [Reference](boundary-conditions.md) |
| `tensorDivDivMatrix` | [Reference](poisson.md) |
| `tensorField` | [Reference](fields.md) |
| `tensorLaplaceMatrix` | [Reference](poisson.md) |
| `tensorToMatrix` | [Reference](constitutive.md) |
| `tensorToVector` | [Reference](constitutive.md) |
| `thermalLoadVector` | [Reference](heat.md) |
| `trace` | [Reference](operators.md) |
| `traceLaplaceMatrix` | [Reference](poisson.md) |
| `unitTensor` | [Reference](operators.md) |
| `updateContact!` | [Reference](contact.md) |
| `vectorField` | [Reference](fields.md) |
| `ε` | [Reference](multifield.md) |
| `∂t` | [Reference](field-calculus.md) |
| `∂x` | [Reference](field-calculus.md) |
| `∂y` | [Reference](field-calculus.md) |
| `∂z` | [Reference](field-calculus.md) |
| `∇` | [Reference](field-calculus.md) |
| `∇_old` | [Reference](nonlinear.md) |
| `∫` | [Reference](multifield.md) |
| `∫N_c_N_dΩ` | [Reference](poisson.md) |
| `∫N_c_dΩ` | [Reference](linear.md) |
| `∫N_c_∂N∂x_dΩ` | [Reference](poisson.md) |
| `∫N_c_∂N∂y_dΩ` | [Reference](poisson.md) |
| `∫N_c_∂N∂z_dΩ` | [Reference](poisson.md) |
| `∫Γ` | [Reference](multifield.md) |
| `∫Ω` | [Reference](multifield.md) |
| `∫∇N_c_∇N_dΩ` | [Reference](poisson.md) |
| `∫∇oN_c_∇oN_dΩ` | [Reference](poisson.md) |
| `∫∇xN_c_∇xN_dΩ` | [Reference](poisson.md) |
| `⋅` | [Reference](operators.md) |
