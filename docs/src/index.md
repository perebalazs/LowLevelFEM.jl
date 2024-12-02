# LowLevelFEM

Solution of a problem in linear elasticity using Finite Element Method consists of solution of the stiffness matrix **K** and load vector **f**, modifying them according to the boundary conditions (getting **K**' and **f**'), solving the displacement field **q**' as the result of the system of equations **K**'**q**=**f**', solving the stress field from **q** and visualize them.
The above described steps can be easily performed using the LowLevelFEM package. Each step means a function with the appropriate parameters, while at any step it is possible to perform an arbitrary operation with the quantities calculated in the meantime. For example the strain energy can be solved as U=1/2**q**^T^**Kq**, for which the code is simply ```U=q'*K*q/2.```(see Examples)

## Features

- Sketching the geometry, making the FE mesh with [GMSH](https://gmsh.info).
- Solving problems from linear elasticity,

  - 2D problems,
    - Plane stress,
    - Plane strain,
    - Axisymmetric,
  - 3D problem (solid body),
- which means the creation of the stiffness matrix **K** and the mass matrix **M** of the problem using arbitrary

  - element types (line, triangle, rectangle, tetrahedron, hexahedron, pyramid, wedge)
  - approximation order (up to ten, Lagrange polynomials)
- Applying

  - distributed forces on arbitrary *physical groups* (see [GMSH](https://gmsh.info)),

    - Lines (in 2D: surface force, in 3D: edge force),
    - Surfaces (in 2D: body force, in 3D: traction),
    - Volumes (in 3D: body force),
  - concentrated forces on nodes,

    which means the calculation of the load vector **f**.
- Constraints on physical groups (nodes on points, edges, surfaces and volumes).
- Elastic support
- Giving loads as functions
- Giving displacement constraints as functions
- Different materials on each physical group
- Solves stress, stain and heat flux field as element result (possibly jumps at the element boundaries) or as nodal results.
- Resultant of "load vector" type quantities on arbitrary physical group (in GMSH).
- Applying initial conditions (displacement and velocity) on arbitrary points, edges, surfaces, volumes and on combinations of them.
- Solution of static and dynamic (transient with central difference method, Newmark and HHT-α) problems,
- Displaying the results (scalar or vector displacements, scalar or tensor stresses and strains) with [GMSH](https://gmsh.info).
  - When dynamic problems are solved animations are also possible (click on $\triangleright$).
- Plotting arbitrary results on paths.
- Solves the damping matrix of structures in case of proportional damping
  - using Rayleigh-damping (**C**=α**M**+β**K**) or
  - using Caughey-damping (**C**=α**M**+β₁**K**+β₂**KM⁻¹K**+β₃**KM⁻¹KM⁻¹K**+⋅⋅⋅).
- Solves the stability analysis transient problems (spectral radius, period error, physical damping ratio, algorithmic damping ratio)
- Heat conduction problems
    - Conductivity and heat capacity matrices,
    - Temperature boundary conditions (also with functions)
    - Loads:
        * Heat flux on boundaries (also with functions)
        * Heat source inside the bodies (also with functions)
        * Heat convection
    - Heat expansion
    - Thermal loading in stress analysis (thermal stresses)
- Modal analysis (eigenfrequencies, modal shapes).

## Planned features

- [ ] 3D (and  2D) truss structures
- [ ] 3D (and 2D) beam structures
- [ ] Shells
- [ ] MultiPoint Constraint (like MPC184 in Ansys)
- [ ] Contact problems,

  - [ ] in 2D,
  - [ ] in 3D,
  - [ ] with penalty method
  - [ ] with Lagrange multiplier method.
- [ ] Defining and using coordinate systems,

  - [ ] cartesian at arbitrary position and arbitrary orientation,
  - [ ] cylindrical.
- [ ] Defining displacement initial condition as a function of x, y and z.
- [ ] Defining velocity initial condition as a function of x, y and z.
- [ ] Transient problems in heat conduction.
- [ ] Finite deformations.
- [ ] Linear buckling.

Any suggestions are welcome.
