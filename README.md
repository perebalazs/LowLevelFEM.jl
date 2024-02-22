[![Build Status](https://github.com/perebalazs/LowLevelFEM.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/perebalazs/LowLevelFEM.jl/actions/workflows/CI.yml?query=branch%3Amain)

# LowLevelFEM

LowLevelFEM package gives the possibility to solve finite element problems using matrix equations. For example to solve problems from linear elasticity the stiffness matrix $\mathbf{K}$ and load vector $\mathbf{f}$ are needed. After applying the neccesery boundary conditions to $\mathbf{K}$ and $\mathbf{f}$ (and getting $\tilde{\mathbf{K}}$ and $\tilde{\mathbf{f}}$) the system of linear equations $\tilde{\mathbf{K}}\mathbf{q}=\tilde{\mathbf{f}}$ have to be solved. With the functions of the package the above described tasks can be performed easily.

# Features

- Sketching the geometry, making the FE mesh with [GMSH](https://gmsh.info).
- Solving problems from linear elasticity,
  - 2D problems,
    - Plane stress,
    - Plane strain,
  - 3D problem (solid body),
- which means the creation of the stiffness matrix $\mathbf{K}$ of the problem using arbitrary
  - element types (line, triangle, rectangle, tetrahedron, hexahedron, pyramid, wedge)
  - approximation order (up to ten, Lagrange polynomials)
- Applying distributed forces on arbitrary *physical groups* (see [GMSH](https://gmsh.info)),
  - Lines (in 2D: surface force, in 3D: edge force),
  - Surfaces (in 2D: body force, in 3D: traction),
  - Volumes (in 3D: body force),
- concentrated forces on nodes,
- which means the calculation of the load vector $\mathbf{f}$.
- Constraints on physical groups (nodes on points, edges, surfaces and volumes).
- Applying initial conditions on arbitrary points, edges, surfaces, volumes and on combinations of them.
- Solution of static and dynamic (transient with central difference method) problems,
- which means the generations of the mass matrix $\mathbf{M}$.
- Displaying the results (scalar or vector displacements, scalar or tensor stresses) with [GMSH](https://gmsh.info).
  - When dynamic problems are solved animations are also possible (click on $\triangleright$).
- Plotting arbitrary results on paths.

# Planed features

- [ ] 2D axisymmetric problems
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

# Examples

## 2D Cantilever

cantilever2D.jl
```Julia
using LinearAlgebra, SparseArrays
import LowLevelFEM as FEM
using LowLevelFEM

gmsh.initialize()

gmsh.open("cantilever2D.geo")
problem = FEM.Problem(type="PlaneStress")

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
using LinearAlgebra, SparseArrays
import LowLevelFEM as FEM
using LowLevelFEM

gmsh.initialize()

gmsh.open("cantilever3D.geo")
problem = FEM.Problem()

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

---

