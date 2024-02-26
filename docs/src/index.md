# LowLevelFEM

Solution of a problem in linear elasticity using Finite Element Method consists of solution of the stiffness matrix $\mathbf{K}$ and load vector $\mathbf{f}$, modifying them according to the boundary conditions (getting $\tilde{\mathbf{K}}$ and $\tilde{\mathbf{f}}$), solving the displacement field $\mathbf{q}$ as the result of system of equations $\tilde{\mathbf{K}}\mathbf{q}=\tilde{\mathbf{f}}$, solving the stress field from $\mathbf{q}$ and visualize them.
The above described steps can be easily performed using the LowLevelFEM package. Each step means a function with the appropriate parameters, while at any step it is possible to perform an arbitrary operation with the quantities calculated in the meantime. For example the strain energy can be solved as $U=\frac{1}{2}\mathbf{q}^T\mathbf{K}\mathbf{q}$, for which the code is simply ```U=q'*K*q/2.```(see Examples)

## Features

- Sketching the geometry, making the FE mesh with [GMSH](https://gmsh.info).
- Solving problems from linear elasticity,
  
  - 2D problems,
    - Plane stress,
    - Plane strain,
  - 3D problem (solid body),
- which means the creation of the stiffness matrix $\mathbf{K}$ of the problem using arbitrary
  
  - element types (line, triangle, rectangle, tetrahedron, hexahedron, pyramid, wedge)
  - approximation order (up to ten, Lagrange polynomials)
- Applying
  
  - distributed forces on arbitrary *physical groups* (see [GMSH](https://gmsh.info)),
    
    - Lines (in 2D: surface force, in 3D: edge force),
    - Surfaces (in 2D: body force, in 3D: traction),
    - Volumes (in 3D: body force),
  - concentrated forces on nodes,
    which means the calculation of the load vector $\mathbf{f}$.
- Constraints on physical groups (nodes on points, edges, surfaces and volumes).
- Applying initial conditions on arbitrary points, edges, surfaces, volumes and on combinations of them.
- Solution of static and dynamic (transient with central difference method) problems,
- which means the generations of the mass matrix $\mathbf{M}$.
- Displaying the results (scalar or vector displacements, scalar or tensor stresses) with [GMSH](https://gmsh.info).
  
  - When dynamic problems are solved animations are also possible (click on $\triangleright$).
- Plotting arbitrary results on paths.

## Planned features

- [ ] 2D axisymmetric problem
- [ ] 3D (and 2D) beam structures
- [ ] Shells
- [ ] Giving loads and prescribed displacements with functions
- [ ] Different material properties on physical groups
- [ ] Contact problems,
  
  - [ ] in 2D,
  - [ ] in 3D,
  - [ ] with Lagrange multiplier method.
- [ ] Using different materials within a model.
- [ ] Defining and using coordinate systems,
  
  - [ ] cartesian at arbitrary position and arbitrary orientation,
  - [ ] cylindrical.
- [ ] Finite deformations.
- [ ] Heat conduction problems,
  
  - [ ] solving conductivity matrix,
  - [ ] solving heat capacity matrix.
- [ ] Dynamic transient problems with HHT-α (or Newmark).
- [ ] Linear buckling.
- [ ] Modal analysis (eigenfrequencies, modal shapes).

Any suggestions are welcome.

