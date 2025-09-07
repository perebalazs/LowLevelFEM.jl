
---
title: 'LowLevelFEM.jl: A lightweight finite element toolbox in Julia'
tags:
  - Julia
  - finite element method
  - structural mechanics
  - continuum mechanics
  - open source
authors:
  - name: Balázs Pere
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: Department of Applied Mechanics, Széchenyi István University, Győr, Hungary
    index: 1
date: 2025-09-07
bibliography: paper.bib
---

# Summary

LowLevelFEM.jl is a finite element method (FEM) toolbox written entirely in the Julia
programming language. Its design philosophy emphasizes **simplicity, transparency, and
performance**, making it suitable for both educational and research purposes in mechanics
and engineering. Unlike large frameworks that rely on domain-specific languages or
compiled backends, LowLevelFEM provides users with direct access to FEM building blocks
in Julia, enabling full control over discretization, assembly, and solution steps.

The package currently supports two- and three-dimensional solid mechanics problems
including plane stress, plane strain, axisymmetric, and 3D solid analyses. Its minimalistic
design lowers the entry barrier for students while still offering enough flexibility for
advanced users to **inspect and control algorithms step by step, process intermediate
results during the solution procedure, and perform operations with scalar, vector, and
tensor fields**. LowLevelFEM uses **Gmsh** as its pre- and post-processor, ensuring
compatibility with a widely adopted meshing and visualization tool.

Thanks to Julia’s just-in-time compilation and multiple dispatch, the code remains concise
while achieving performance comparable to traditional FEM codes in C/C++ or Fortran.
LowLevelFEM is released under the MIT license and distributed via the Julia General
registry. Documentation, tutorials, and examples are available at
<https://perebalazs.github.io/LowLevelFEM.jl/stable/>.

## Example

Below is a simple example illustrating a typical LowLevelFEM workflow using Gmsh for
pre- and post-processing:

```julia
using LowLevelFEM

# `gmsh` is exported by LowLevelFEM
gmsh.initialize()
gmsh.open("your_model.geo")

mat = material("body", E=2e5, ν=0.3)
prob = Problem([mat], type=:PlaneStress)  # :Solid, :PlaneStrain, :AxiSymmetric, :HeatConduction, ...

bc   = displacementConstraint("supp", ux=0, uy=0)
force = load("load", fy=-1)

q = solveDisplacement(prob, [force], [bc])
S = solveStress(q)

showDoFResults(q, :uvec)
showStressResults(S, :s)

openPostProcessor()
gmsh.finalize()
````

Alternatively, a lower-level sequence:

```julia
K = stiffnessMatrix(prob)
f = loadVector(prob, [force])
applyBoundaryConditions!(K, f, [bc])
q = K \ f

# Simple Hooke's law stress computation
A = (u ∘ ∇ + ∇ ∘ u) / 2
I = unitTensor(A)
S = E / (1 + ν) * (A + ν / (1 - 2 ν) * trace(A) * I)
```

# Statement of need

Finite element simulations are essential in many fields of engineering, especially in
solid mechanics and structural analysis. However, educational and research communities
often face two challenges:

1. **Accessibility**: Commercial FEM packages are expensive and closed-source,
   limiting their use in academic teaching and reproducible research.
2. **Extensibility**: Large open-source frameworks such as FEniCS or deal.II provide
   powerful high-level interfaces but are difficult to extend at the low-level assembly
   stage without diving into C++ backends.

LowLevelFEM addresses these challenges by offering a **lightweight Julia-only
implementation** that exposes all the core FEM routines directly in the high-level
language. This makes the package particularly well-suited for:

* Teaching FEM concepts in undergraduate and graduate courses.
* Rapid prototyping of new FEM formulations and **non-standard algorithms**.
* Research projects where step-by-step inspection of the solution process and manipulation
  of intermediate fields is required.
* Demonstrations of Julia’s potential as a performant and expressive language for
  numerical mechanics.

By combining transparent algorithms with Julia’s scientific ecosystem (e.g.
LinearAlgebra.jl for scalar, vector and tensor field operations, Plots.jl for visualization), and by relying
on **Gmsh** for pre- and post-processing, LowLevelFEM serves as a bridge between pedagogy
and advanced research workflows.

# References

* Bezanson, J., Edelman, A., Karpinski, S., & Shah, V. B. (2017). Julia: A fresh approach to numerical computing. *SIAM Review*, 59(1), 65–98. [https://doi.org/10.1137/141000671](https://doi.org/10.1137/141000671)

* Alnæs, M., Blechta, J., Hake, J., Johansson, A., Kehlet, B., Logg, A., Richardson, C., et al. (2015). The FEniCS Project Version 1.5. *Archive of Numerical Software*, 3(100), 9–23. [https://doi.org/10.11588/ans.2015.100.20553](https://doi.org/10.11588/ans.2015.100.20553)

* Bangerth, W., Hartmann, R., & Kanschat, G. (2007). deal.II – A general-purpose object-oriented finite element library. *ACM Transactions on Mathematical Software*, 33(4). [https://doi.org/10.1145/1268776.1268779](https://doi.org/10.1145/1268776.1268779)

* Badia, S., & Verdugo, F. (2020). Gridap: An extensible Finite Element toolbox in Julia. *Journal of Open Source Software*, 5(52), 2520. [https://doi.org/10.21105/joss.02520](https://doi.org/10.21105/joss.02520)

* Geuzaine, C., & Remacle, J.-F. (2009). Gmsh: A three-dimensional finite element mesh generator with built-in pre- and post-processing facilities. *International Journal for Numerical Methods in Engineering*, 79(11), 1309–1331. [https://doi.org/10.1002/nme.2579](https://doi.org/10.1002/nme.2579)