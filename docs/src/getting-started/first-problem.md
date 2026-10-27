# First Problem

Solve a 3D elastic cube with weak-form assembly, then compute strain and stress
with field calculus. The tiny mesh keeps the example quick; engineering accuracy
requires a convergence study for the actual geometry and loading.

```julia
# docs-check: first-elasticity
using LowLevelFEM, LinearAlgebra

gmsh.initialize()
try
    structured_box_mesh(n=1, order=1)
    mat = Material("body", E=2e5, ν=0.3)
    U = Field([mat], type=:VectorField, dim=3, fieldName=:u)
    K = ∫(SymGrad(U) ⋅ D(:Solid, mat) ⋅ SymGrad(U))
    f = ∫(U ⋅ [1.0, 0.0, 0.0]; Γ="right")
    bc = BoundaryCondition("left", field=U, ux=0, uy=0, uz=0)
    u = solveField(Symmetric(K), f; support=[bc])

    εtensor = (u ∘ ∇ + ∇ ∘ u) / 2
    σ = mat.E/(1+mat.ν) *
        (εtensor + mat.ν/(1-2mat.ν)*trace(εtensor)*unitTensor(εtensor))
    @assert all(isfinite, u.a)
    @assert all(A -> all(isfinite, A), σ.A)
    showDoFResults(u)
    showStressResults(σ)
    # openPostProcessor()  # optional interactive GUI
finally
    gmsh.finalize()
end
```

The helper creates `"body"`, `"left"` and `"right"` physical groups.
`SymGrad(U)` is a DSL operator; `∇(u)` differentiates the solved field. Natural
loads enter the RHS integral and essential conditions enter the solver.
To inspect Gmsh views interactively, open the postprocessor before finalizing
the session.

## Next steps

- [Scalar diffusion with an exact solution](../manual/assembly-and-solvers.md)
- [Mesh and physical groups](mesh-and-physical-groups.md)
- [Fields and results](fields-and-results.md)
- [Tutorials](../tutorials/index.md)
