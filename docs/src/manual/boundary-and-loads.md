# Boundary Conditions and Loads

Boundary conditions and loads are defined independently of the finite element assembly. They can be created from constants, functions, or physical groups defined in the Gmsh model, then passed directly to the solver.

A typical workflow is:

```julia
bc = displacementConstraint("left", ux=0, uy=0)

force = load("right", fy=-1)

u = solveDisplacement(problem,
    support=[bc],
    load=[force])
```

Supported boundary conditions and loads include:

- displacement and temperature constraints
- nodal and distributed mechanical loads
- body forces and heat sources
- heat flux and convection boundary conditions
- elastic supports
- initial displacement, velocity and temperature fields

For detailed function descriptions and additional options, see the API reference below.

## Related API

- [Reference: Preprocessing](../reference/preprocessing.md)
- [Reference: Heat](../reference/heat.md)
