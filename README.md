[![Build Status](https://github.com/perebalazs/LowLevelFEM.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/perebalazs/LowLevelFEM.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://perebalazs.github.io/LowLevelFEM.jl/dev)

# LowLevelFEM

Solution of a problem in linear elasticity using Finite Element Method consists of solution of the stiffness matrix $\mathbf{K}$ and load vector $\mathbf{f}$, modifying them according to the boundary conditions (getting $\tilde{\mathbf{K}}$ and $\tilde{\mathbf{f}}$), solving the displacement field $\mathbf{q}$ as the result of the system of equations $\tilde{\mathbf{K}}\mathbf{q}=\tilde{\mathbf{f}}$, solving the stress field from $\mathbf{q}$ and visualize them.
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
- Displaying the results (scalar or vector displacements, scalar or tensor stresses and strains) with [GMSH](https://gmsh.info).

  - When dynamic problems are solved animations are also possible (click on $\triangleright$).
- Plotting arbitrary results on paths.
- Solves the damping matrix of structures in case of proportional damping

  - using Rayleigh-damping (**C**=α**M**+β**K**) or
  - using Caughey-damping (**C**=α**M**+β₁**K**+β₂**KM⁻¹K**+β₃**KM⁻¹KM⁻¹K**+⋅⋅⋅).

## Planned features

- [ ] 3D (and  2D) truss structures
- [x] 2D axisymmetric problem
- [ ] 3D (and 2D) beam structures
- [ ] Shells
- [x] Giving loads with functions
- [x] Giving prescribed displacements with functions
- [ ] MultiPoint Constraint (like MPC184 in Ansys)
- [X] Different material properties on physical groups
- [ ] Contact problems,

  - [ ] in 2D,
  - [ ] in 3D,
  - [ ] with Lagrange multiplier method.
- [ ] Defining and using coordinate systems,

  - [ ] cartesian at arbitrary position and arbitrary orientation,
  - [ ] cylindrical.
- [x] Defining load vector as a function of x, y and z.
- [x] Defining displacement boundary conditions as a function of x, y and z.
- [ ] Defining displacement initial condition as a function of x, y and z.
- [ ] Defining velocity initial condition as a function of x, y and z.
- [ ] Finite deformations.
- [ ] Heat conduction problems,

  - [ ] solving conductivity matrix,
  - [ ] solving heat capacity matrix.
- [x] Dynamic transient problems with HHT-α (or Newmark).
- [ ] Linear buckling.
- [ ] Modal analysis (eigenfrequencies, modal shapes).

Any suggestions are welcome.

# Examples

## 2D Cantilever

cantilever2D.jl

```Julia
import LowLevelFEM as FEM
using LowLevelFEM

gmsh.initialize()

gmsh.open("cantilever2D.geo")
mat = FEM.material("body", E=2.e5, ν=0.3)
problem = FEM.Problem([mat], type="PlaneStress")

supp = FEM.displacementConstraint("supp", ux=0, uy=0)
load = FEM.load("load", fy=-1)

K = FEM.stiffnessMatrix(problem)
f = FEM.loadVector(problem, [load])

FEM.applyBoundaryConditions!(problem, K, f, [supp])

q = FEM.solveDisplacement(K, f)
S = FEM.solveStress(problem, q)

u = FEM.showDoFResults(problem, q, "uvec", name="uvec", visible=false)
ux = FEM.showDoFResults(problem, q, "ux", name="ux", visible=false)
uy = FEM.showDoFResults(problem, q, "uy", name="uy", visible=false)

s = FEM.showStressResults(problem, S, "s", name="σ", visible=true, smooth=true)
sx = FEM.showStressResults(problem, S, "sx", name="σx", visible=false, smooth=true)
sy = FEM.showStressResults(problem, S, "sy", name="σy", visible=false, smooth=true)
sxy = FEM.showStressResults(problem, S, "sxy", name="τxy", visible=false, smooth=true)

FEM.plotOnPath(problem, "path", sx, 100, name="σx", visible=false);
FEM.plotOnPath(problem, "path", sxy, 100, name="τxy", visible=false);
FEM.plotOnPath(problem, "path", ux, 100, name="ux", visible=false);

gmsh.fltk.run()
gmsh.finalize()
```

cantilever2D.geo

```gmsh
SetFactory("OpenCASCADE");

Rectangle(1) = {0, 0, 0, 100, 10, 0};

Physical Curve("supp", 5) = {4};
Physical Curve("load", 6) = {2};
Physical Surface("body", 7) = {1};

Recombine Surface {1};

Transfinite Line {2,4} = 4;
Transfinite Line {1,3} = 31;
Transfinite Surface {1};

Mesh.ElementOrder = 3;

SetName "cantilever2D";
Mesh 2;

Point(5) = {10, 0, 0, 1.0};
Point(6) = {10, 10, 0, 1.0};
Line(5) = {5, 6};

Physical Curve("path", 8) = {5};
```

## 3D Cantilever

cantilever3D.jl

```Julia
import LowLevelFEM as FEM
using LowLevelFEM

gmsh.initialize()

gmsh.open("cantilever3D.geo")
mat = FEM.material("body", E=2.e5, ν=0.3)
problem = FEM.Problem([mat])

supp = FEM.displacementConstraint("supp", ux=0, uy=0, uz=0)
load = FEM.load("load", fy=-1)

K = FEM.stiffnessMatrix(problem)
f = FEM.loadVector(problem, [load])

FEM.applyBoundaryConditions!(problem, K, f, [supp])

q = FEM.solveDisplacement(K, f)
S = FEM.solveStress(problem, q)

u = FEM.showDoFResults(problem, q, "uvec", name="uvec", visible=false)
ux = FEM.showDoFResults(problem, q, "ux", name="ux", visible=false)
uy = FEM.showDoFResults(problem, q, "uy", name="uy", visible=false)
uz = FEM.showDoFResults(problem, q, "uz", name="uz", visible=false)

s = FEM.showStressResults(problem, S, "s", name="σ", visible=true, smooth=true)
sx = FEM.showStressResults(problem, S, "sx", name="σx", visible=false, smooth=true)
sy = FEM.showStressResults(problem, S, "sy", name="σy", visible=false, smooth=true)
sz = FEM.showStressResults(problem, S, "sz", name="σz", visible=false, smooth=true)
sxy = FEM.showStressResults(problem, S, "sxy", name="τxy", visible=false, smooth=true)
syz = FEM.showStressResults(problem, S, "syz", name="τyz", visible=false, smooth=true)
szx = FEM.showStressResults(problem, S, "szx", name="τzx", visible=false, smooth=true)

FEM.plotOnPath(problem, "path", sx, 100, name="σx", visible=false);
FEM.plotOnPath(problem, "path", sxy, 100, name="τxy", visible=false);
FEM.plotOnPath(problem, "path", ux, 100, name="ux", visible=false);

gmsh.fltk.run()
gmsh.finalize()
```

cantilever3D.geo

```gmsh
SetFactory("OpenCASCADE");

Box(1) = {0, 0, 0, 100, 10, 10};

Physical Surface("supp", 13) = {1};
Physical Surface("load", 14) = {2};
Physical Volume("body", 15) = {1};

Recombine Surface {1:6};

Transfinite Line {1:8} = 4;
Transfinite Line {9:12} = 31;
Transfinite Surface {1:6};
Transfinite Volume {1};

Mesh.ElementOrder = 3;

SetName "cantilever3D";
Mesh 3;

Point(9) = {10, 0, 5, 1.0};
Point(10) = {10, 10, 5, 1.0};
Line(13) = {9, 10};

Physical Curve("path", 16) = {13};
```

For more examples see [Documentation](https://docs.juliahub.com/General/LowLevelFEM) or [examples on GitHub](https://github.com/perebalazs/LowLevelFEM.jl/tree/main/examples)

---

