export Problem, Field, Material, getEigenValues, material
export displacementConstraint, load, elasticSupport, BoundaryCondition, BoundaryConditionFields
export LoadCondition, MultiPointConstraint, MPC
export temperatureConstraint, heatFlux, heatSource, heatConvection
export field, scalarField, vectorField, tensorField, ScalarField, VectorField, TensorField
export mergeFields
export SystemMatrix, SystemVector, subSystemMatrix
export constrainedDoFs, freeDoFs, allDoFs, DoFs
export elementsToNodes, nodesToElements, elementsToElements
export projectTo2D, expandTo3D, isNodal, isElementwise
export fieldError, resultant, integrate, ∫
export normalVector, normalVectors, tangentVector, tangentVectors
export rotateNodes, CoordinateSystem, Transformation
export showDoFResults, showModalResults, showBucklingResults
export showStrainResults, showStressResults, showElementResults, showHeatFluxResults
export plotOnPath, showOnSurface, plotOnBeam
export openPreProcessor, openPostProcessor, setParameter, setParameters, getParameter
export probe
export saveField, loadField, isSaved
export ∂x, ∂y, ∂z, ∂t
export structured_rect_mesh, structured_box_mesh, line_mesh, openGeometry
export renumberNodes!, ndofs
export ConstitutiveMatrix, D

# Definitions remain in their original evaluation order.
include("general/materials_and_problems.jl")
include("general/system_matrices.jl")
include("general/fields.jl")
include("general/boundary_conditions.jl")
include("general/field_calculus.jl")
include("general/field_builders.jl")
include("general/field_mapping.jl")
include("general/integration.jl")
include("general/postprocessing.jl")
include("general/io.jl")
include("general/mesh_helpers.jl")
