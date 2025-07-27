# LowLevelFEM.jl

**Repository:** https://github.com/perebalazs/LowLevelFEM.jl  
**Purpose:** Solve linear-elasticity finite element problems with fine-grained, low-level control over mesh, assembly, boundary conditions, and visualization. It uses engineering approach, no need to know what is the *weak form*.

---

## 🎯 Key Features

- **Geometry & Mesh Generation** via GMSH (2D and 3D)  
- Assemble **stiffness matrix (K)** and **load vector (f)** for linear elasticity problems  
- Supports arbitrary **element types**: line, triangle, rectangle, tetrahedron, hexahedron, pyramid, wedge  
- **High-order Lagrange basis functions**, up to 10th degree 
- **Heterogeneous materials** by physical groups  
- **Apply distributed and nodal forces** on physical groups (points, lines, surfaces, volumes). Distributes forces can be given as functions of the position.
- **Boundary conditions** on points, curves, surfaces, or volumes  
- **Static and dynamic (transient)** solvers (mass matrix `M`, central-difference or HHT time integration)
- **Modal analysis** to solve eigenfrequencies and modal shapes, even if the srtucture is prestressed.
- **Buckling** analysis for linear elasticity
- **Heat conduction** problems, thermal loads can be heat flux or convection (radiation not yet implemented). Thermal stresses and transient problems also can be computed.
- **Arbitrary coordinale systems** given with user defined function (eg. curvilinear)
- Possibility to give **user defined scalar, vector or tensor fields** as a finite element result in order to compare results with analytic ones.
- **Visualization via GMSH**, including displacement fields, stress fields, and animations for dynamic solutions   
- Plot results along **user-defined paths**

---

## 🎖 Strengths

| Feature               | Description |
|-----------------------|-------------|
| **Fine-grained control** | Users can intervene at every step (mesh → K/f assembly → BC → solve → postprocess), ideal for research and custom workflows |
| **Modular and extensible** | Each phase is a separate function, making it easy to customize or augment assembly or solution steps  |
| **High-order elements** | Supports Lagrange approximations up to 10th order, enabling high-accuracy solutions in both 2D and 3D settings  |
| **Well-documented** | Docstrings with Documenter.jl, automatic docs deployed on GitHub Pages and JuliaHub  |

---

## 🛠 Roadmap

- 2D and 3D **axisymmetric analysis**  
- **Beam and shell elements** (e.g., beams, plates, shells)  
- **Multi-point constraints** (e.g., MPC184-style constraints)  
- **Contact mechanics**, using Lagrange multipliers and penalty method 
- **Finite deformations and nonlinear material laws** eg. neo-Hooke, Mooney-Rivlin, etc. (practically finished, will be released in the next version this summer)

---

## 🔧 Example Usage

```julia
using LowLevelFEM
import LowLevelFEM as FEM

gmsh.initialize()
gmsh.open("cantilever2D.geo")

mat = FEM.material("body", E=2e5, ν=0.3)
prob = FEM.Problem([mat], type=:PlaneStress)

bc = FEM.displacementConstraint("supp", ux=0, uy=0)
load = FEM.load("load", fy=-1)

K = FEM.stiffnessMatrix(prob)
f = FEM.loadVector(prob, [load])

FEM.applyBoundaryConditions!(prob, K, f, [bc])

q = FEM.solveDisplacement(K, f)
S = FEM.solveStress(prob, q)

u = FEM.showDoFResults(prob, q, :uvec)
s = showStressResults(prob, S, :s)

plotOnPath(prob, "path", s)
FEM.openPostProcessor()
gmsh.finalize()
```

This workflow sketch demonstrates geometry loading, problem setup, assembly, solver execution, postprocessing, and visualization. For the sake of shorter scripts a few rows can be combined: `q = FEM.solveDisplacement(prob, [load], [supp])`

---

## 🔗 Ecosystem & Comparisons

- **juliaPDE** – https://github.com/JuliaPDE/SurveyofPDEPackages?tab=readme-ov-file#fem

---

## ℹ️ Summary

**LowLevelFEM.jl** is a powerful, flexible, and modular Julia package for solving linear-elasticity FEM problems with full control at every step. Ideal for researchers, educators, and practitioners who require customizable assembly, solver, and visualization workflows.

---

