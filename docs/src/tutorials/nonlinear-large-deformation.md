# Nonlinear Formulations with the DSL

A nonlinear model needs kinematics, a constitutive response, an assembled
residual/tangent, and an equilibrium iteration. The DSL exposes the assembly
pieces; `solveField` solves a linear system and does not run Newton iteration.
The old dedicated total-Lagrangian assembly functions remain in the
[legacy nonlinear reference](../reference/nonlinear.md).

## Kinematics from a displacement field

For a three-component displacement field `u`, construct the deformation gradient
and right Cauchy–Green tensor:

```julia
H = u ∘ ∇
F = unitTensor(H) + H
C = F' * F
```

This is field calculus on computed values. A 2D displacement needs a deliberate
3D embedding and an appropriate plane-strain or plane-stress constitutive
assumption; expansion alone does not solve a plane-stress reduction.

## Energy-based response

Load `Tensors` and define a differentiable `ψ(C,p)` using its tensor operations.
Then evaluate stress and material tangent at element nodes:

```julia
response = energyMaterial(ψ, C; params=(μ=mat.μ, λ=mat.λ))
S = response.S
D_raw = response.D
pvoigt = [1, 2, 3, 6, 4, 5]
D_dsl = D_raw[pvoigt, pvoigt]
```

The permutation matches the current extension's tangent order to the DSL Voigt
order. Use `IIPiolaKirchhoff` or `tangentMatrix` when only one response is needed.
Constitutive sampling and Gauss-point interpolation are distinct operations;
convert coefficient representations where required by the chosen assembly path.
See [Constitutive API](../reference/constitutive.md) for conventions and restrictions.

## Equilibrium and tangent consistency

At each load step, evaluate the residual and a consistent tangent for the chosen
reference/current configuration. Material and geometric terms, follower loads,
constraints and any contact active-set updates must all match that formulation.
Solve the linearized correction with `solveField`, update the displacement,
and check both residual and increment convergence. This is an algorithmic
outline, not a complete generic nonlinear solver.

For a concrete DSL geometric-stiffness application, see
[Linear Buckling and Prestressed Natural Frequencies](buckling-and-prestressed-modes.md).
That tutorial linearizes around a preloaded state; it does not replace a
finite-deformation equilibrium solve.

## Check anisotropic shear ordering independently

This optional material-response example checks the documented tangent
permutation against three distinct engineering shear moduli. It is a
constitutive check, not a nonlinear equilibrium simulation. `Tensors` must
already be installed and is loaded explicitly.

```julia
# docs-check: energy-material; optional
using LowLevelFEM, Tensors

gmsh.initialize()
try
    structured_box_mesh(n=1, order=1)
    P = Field([Material("body")], type=:TensorField, dim=3)
    C = TensorField(P, "body", [1.0 0.2 0.3; 0.2 1.0 0.4; 0.3 0.4 1.0])
    psi(C,p) = p.a/8*(C[1,2]+C[2,1])^2 +
               p.b/8*(C[2,3]+C[3,2])^2 + p.c/8*(C[3,1]+C[1,3])^2
    response = energyMaterial(psi, C; params=(a=2.0,b=3.0,c=5.0))
    perm = [1,2,3,6,4,5]
    D_dsl = response.D[perm,perm]
    for (i,mu) in zip(4:6,(2.0,3.0,5.0))
        @assert all(A -> maximum(abs.(A .- mu)) < 1e-10, D_dsl[i,i].A)
    end
    @assert all(A -> maximum(abs.(A .- 0.4)) < 1e-10, response.S[1,2].A)
finally
    gmsh.finalize()
end
```
