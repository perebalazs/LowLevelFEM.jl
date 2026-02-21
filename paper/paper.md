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
    orcid: 0000-0002-1161-4206
    affiliation: 1
affiliations:
  - name: Department of Applied Mechanics, Széchenyi István University, Győr, Hungary
    index: 1
date: 2025-09-07
bibliography: paper.bib
---

# Summary

LowLevelFEM.jl is a finite element method (FEM) [@zienkiewicz2005] toolbox written entirely in the Julia programming language [@bezanson2017]. Its design philosophy emphasizes **simplicity, transparency, and performance**, making it suitable for both educational and research purposes in mechanics and engineering. Unlike large frameworks that rely on domain-specific languages or compiled backends, LowLevelFEM provides users with direct access to FEM building blocks in Julia, enabling full control over discretization, assembly, and solution steps.

The package currently supports two- and three-dimensional solid mechanics problems including plane stress, plane strain, axisymmetric, and 3D solid analyses. Its minimalistic design lowers the entry barrier for students while still offering enough flexibility for advanced users to **inspect and control algorithms step by step, process intermediate results during the solution procedure, and perform operations with scalar, vector, and tensor fields**. LowLevelFEM uses **Gmsh** [@geuzaine2009] as its pre- and post-processor, ensuring compatibility with a widely adopted meshing and visualization tool. The toolbox also supports geometrically nonlinear formulations based on a Total Lagrangian framework and energy-driven constitutive modeling, enabling research-level investigations beyond linear elasticity.

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

u = solveDisplacement(prob, load=[force], support=[bc])
S = solveStress(u)

showDoFResults(u)
showDoFResults(u, :ux)
showStressResults(S)
showStressResults(S, :sxy, name="Shear stress")

openPostProcessor()
gmsh.finalize()
```

Note: physical group names in the geometry (created in Gmsh) must match the strings used above (e.g., "body", "supp", "load").

Alternatively, a lower-level sequence:

```julia
K = stiffnessMatrix(prob)
f = loadVector(prob, [force])
u = solveDisplacement(K, f, support=[bc])

# Simple Hooke's law stress computation
E = mat.E
ν = mat.ν
A = (u ∘ ∇ + ∇ ∘ u) / 2
I = TensorField(prob, "body", [1 0 0; 0 1 0; 0 0 1])
S = E / (1 + ν) * (A + ν / (1 - 2ν) * trace(A) * I)
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

# State of the field

Several finite element frameworks are available in the Julia ecosystem and beyond. High-level Julia packages such as Gridap.jl [@badia2020] provide domain-specific abstractions for variational formulations, while Ferrite.jl and JuAFEM.jl focus on flexible but more structured implementations. Outside Julia, frameworks such as FEniCS [@logg2012] and deal.II [@bangerth2015] offer powerful and mature solutions implemented primarily in C++.  

LowLevelFEM differs from these frameworks in its explicit low-level design philosophy. Instead of abstracting away the assembly process through domain-specific languages or symbolic formulations, it exposes stiffness matrix construction, load vector assembly, and boundary condition handling directly in Julia code. This makes it particularly suitable for educational purposes and for researchers who need full algorithmic transparency and control over intermediate computational steps.  

# Software design

LowLevelFEM is structured around a minimal set of core abstractions. The `Problem` type encapsulates material definitions, physical groups imported from Gmsh, and analysis settings (e.g., plane stress, 3D solid, heat conduction). Assembly routines such as `stiffnessMatrix`, `loadVector`, and `applyBoundaryConditions!` operate directly on these problem definitions.  

The package separates:  

- mesh and geometry handling (delegated to Gmsh),  
- operator assembly,  
- solution procedures,  
- post-processing of scalar, vector, and tensor fields.  

This modular structure allows users to either follow a high-level workflow (`solveDisplacement`, `solveStress`) or construct custom pipelines at a lower level. The implementation leverages Julia’s multiple dispatch and just-in-time compilation to maintain readable code while achieving competitive performance.

In addition to linear formulations, the software includes a Total Lagrangian nonlinear pipeline based on strain energy functions. Constitutive models can be defined through free energy densities, from which stress measures and corresponding tangent operators are automatically derived. This operator-oriented design makes it possible to experiment with alternative bilinear and nonlinear forms without modifying compiled backends.  

The package provides algebraic and differential operations on scalar, vector, and tensor fields, allowing users to compose continuum-mechanics expressions directly in Julia syntax, including, for example, gradient, divergence, tensor contraction, and trace operations defined at the discrete level. This feature is particularly valuable for educational demonstrations and rapid prototyping of new formulations.

---

# Research impact statement

LowLevelFEM supports both educational and research activities in computational mechanics. In teaching, it enables students to inspect finite element algorithms step by step without switching languages or interacting with opaque compiled backends. In research, it facilitates rapid prototyping of new constitutive models, nonlinear formulations, and custom operators.  

The package has been used in undergraduate and graduate mechanics courses and serves as a foundation for ongoing developments in nonlinear and multi-field finite element formulations. Its open-source MIT license and integration with the Julia ecosystem promote reproducible research and extensibility. Ongoing developments aim to extend the framework toward multi-field finite element formulations while preserving the same transparent operator-level philosophy.

# AI usage disclosure

Generative AI tools were used for minor language editing and documentation refinement. All scientific concepts, algorithms, and software implementations were designed, developed, and verified by the author.

---

# References
