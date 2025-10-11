---
title: 'LowLevelFEM.jl: A lightweight finite element toolbox in Julia'
tags:

* Julia
* finite element method
* structural mechanics
* continuum mechanics
* open source
  authors:
* name: Balázs Pere
  orcid: 0000-0002-1161-4206
  affiliation: 1
  affiliations:
* name: Department of Applied Mechanics, Széchenyi István University, Győr, Hungary
  index: 1
  date: 2025-09-07
  bibliography: paper.bib

---

# Summary

LowLevelFEM.jl is a finite element method (FEM) [@zienkiewicz2005] toolbox written entirely in the Julia programming language [@bezanson2017]. Its design philosophy emphasizes **simplicity, transparency, and performance**, making it suitable for both educational and research purposes in mechanics and engineering. Unlike large frameworks that rely on domain-specific languages or compiled backends, LowLevelFEM provides users with direct access to FEM building blocks in Julia, enabling full control over discretization, assembly, and solution steps.

The package currently supports two- and three-dimensional solid mechanics problems including plane stress, plane strain, axisymmetric, and 3D solid analyses. Its minimalistic design lowers the entry barrier for students while still offering enough flexibility for advanced users to **inspect and control algorithms step by step, process intermediate results during the solution procedure, and perform operations with scalar, vector, and tensor fields**. LowLevelFEM uses **Gmsh** [@geuzaine2009] as its pre- and post-processor, ensuring compatibility with a widely adopted meshing and visualization tool.

Thanks to Julia’s just-in-time compilation and multiple dispatch, the code remains concise while achieving performance comparable to traditional FEM codes in C/C++ or Fortran. LowLevelFEM is released under the MIT license and distributed via the Julia General registry. Documentation, tutorials, and examples are available at [https://juliahub.com/ui/Packages/General/LowLevelFEM](https://juliahub.com/ui/Packages/General/LowLevelFEM) or [https://perebalazs.github.io/LowLevelFEM.jl/stable/](https://perebalazs.github.io/LowLevelFEM.jl/stable/).

## Example

Below is a simple example illustrating a typical LowLevelFEM workflow using Gmsh for pre- and post-processing:

```julia
using LowLevelFEM

# `gmsh` is exported by LowLevelFEM
gmsh.initialize()
gmsh.open("model.geo")

mat = material("body", E=2e5, ν=0.3)
prob = Problem([mat], type=:PlaneStress)  # :Solid, :PlaneStrain, :AxiSymmetric, :HeatConduction, ...

bc   = displacementConstraint("supp", ux=0, uy=0)
force = load("load", fy=-1)

u = solveDisplacement(prob, [force], [bc])
S = solveStress(u)

showDoFResults(u)
showStressResults(S)

openPostProcessor()
gmsh.finalize()
```

Note: physical group names in the geometry (created in Gmsh) must match the strings used above (e.g., "body", "supp", "load").

Alternatively, a lower-level sequence:

```julia
K = stiffnessMatrix(prob)
f = loadVector(prob, [force])
applyBoundaryConditions!(K, f, [bc])
u = K \ f

# Simple Hooke's law stress computation
A = (u ∘ ∇ + ∇ ∘ u) / 2
I = unitTensor(A)
S = E / (1 + ν) * (A + ν / (1 - 2 ν) * trace(A) * I)
```

# Statement of need

Finite element simulations are essential in many fields of engineering, especially in solid mechanics and structural analysis. However, educational and research communities often face two challenges:

1. **Accessibility**: Commercial FEM packages are expensive and closed-source, limiting their use in academic teaching and reproducible research.
2. **Extensibility**: Large open-source frameworks such as FEniCS [@logg2012] or deal.II [@bangerth2015] provide powerful high-level interfaces but are difficult to extend at the low-level assembly stage without diving into C++ backends.

LowLevelFEM addresses these challenges by offering a **lightweight Julia-only implementation** that exposes all the core FEM routines directly in the high-level language. This makes the package particularly well-suited for:

* Teaching FEM concepts in undergraduate and graduate courses.
* Rapid prototyping of new FEM formulations and **non-standard algorithms**.
* Research projects where step-by-step inspection of the solution process and manipulation of intermediate fields is required.
* Demonstrations of Julia’s potential as a performant and expressive language for numerical mechanics [@bezanson2017].

By combining transparent algorithms with Julia’s scientific ecosystem (e.g. LinearAlgebra.jl, Plots.jl) and by relying on **Gmsh** [@geuzaine2009] for pre- and post-processing, LowLevelFEM serves as a bridge between pedagogy and advanced research workflows. It also complements existing Julia FEM frameworks such as Gridap.jl [@badia2020] and interfaces naturally with linear algebra tools like Arpack.jl [@knyazev2017].

# References
