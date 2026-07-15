# Troubleshooting

This page collects solutions to the most common problems encountered when using LowLevelFEM.

Typical issues include:

- missing or mismatched physical group names
- incompatible problem types and field dimensions
- incorrect boundary-condition or load definitions
- Gmsh initialization or finalization errors

If a simulation does not behave as expected, first verify that:

- all referenced physical groups exist in the Gmsh model,
- the selected problem type matches the mesh and material definition,
- boundary conditions and loads are applied to the intended physical groups.

This section will be expanded as additional common issues are identified.
