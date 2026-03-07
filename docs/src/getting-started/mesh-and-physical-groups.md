# Mesh and Physical Groups

LowLevelFEM identifies model regions through Gmsh physical group names.

## Rule

Strings passed to constructors such as `Material("body", ...)`,
`displacementConstraint("supp", ...)`, and `load("load", ...)`
must match physical group names defined in the mesh.

## Minimal Gmsh Example

```gmsh
Physical Surface("body", 1) = {1};
Physical Curve("supp", 2) = {4};
Physical Curve("load", 3) = {2};
```

## Notes

- Use stable, descriptive names for groups.
- Keep boundary/load groups disjoint where possible for clarity.
- Verify group dimensions (`Ω` vs `Γ`) for weak-form assembly.
