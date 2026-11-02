# Nonlinear Cantilever Dynamics: Secant vs Newton-Type HHT

This tutorial compares two nonlinear iterations for the large-amplitude free
vibration of a cantilever using the weak-form DSL and HHT time integration.

## Notebook

[cantilever_HHT_minimal_secant_vs_Newton.ipynb](https://github.com/perebalazs/LowLevelFEM.jl/blob/main/examples/cantilever_HHT_minimal_secant_vs_Newton.ipynb)

The complete notebook contains both calculations and Gmsh visualization. The
code below follows its cell order and uses only LowLevelFEM.jl.

## Model and initial conditions

The beam is a 2D plane-stress continuum with finite rotations and a
Saint-Venant–Kirchhoff constitutive law based on Green–Lagrange strain. Its
length is `1.0`, height `0.04` and thickness `b = 0.01`, with `E = 1e6`,
`ν = 0.3` and `ρ = 1000.0` in consistent units. The left end is clamped.
The displacement is initially zero; a transverse initial velocity excites the
motion without an external force.

```julia
using LowLevelFEM

structured_rect_mesh(lx=1.0, ly=0.04, y0=-0.02, n=2, order=2)
mat = Material("body"; E=1e6, ν=0.3, ρ=1000.0)
U = Field([mat]; type=:VectorField, dim=2, fieldName=:u)
b, D2 = 0.01, D(:PlaneStress, mat)
support = [BoundaryCondition("left"; field=U, ux=0.0, uy=0.0)]
u0 = initialDisplacement(U, "body"; ux=0.0, uy=0.0)
v0 = initialVelocity(U, "body"; vx=0.0, vy=(x,y,z)->0.25*x^2*(3-x)/2)
M = b*mat.ρ*∫(U⋅U)
K = b*∫(SymGrad(U)⋅D2⋅SymGrad(U))
B0 = [1 0 0 0; 0 0 0 1; 0 1 1 0];
```

## 1. Secant-stiffness iteration

With $F=I+\nabla_Xu$ and $E=(F^TF-I)/2$, the internal force is represented by $f_{\mathrm{int}}(u)=K_s(u)u$. The HHT solver updates $K_s$ at each nonlinear iteration.

```julia
function deformation(state)
    H = ∇(expandTo3D(state.fields[1]))
    F = tensorToMatrix(unitTensor(H) + H)
    B = [F[1,1] 0 F[2,1] 0;
         0 F[1,2] 0 F[2,2];
         F[1,2] F[1,1] F[2,2] F[2,1]]
    return H, B
end

function secant(A, state)
    H, B = deformation(state)
    b*∫(Grad(U)⋅B'⋅D2⋅((B+B0)/2)⋅Grad(U); gauss=4)
end
```

```julia
@time u_sec, v_sec = HHT(K, M, 0*u0, u0, v0, 401, 0.02;
    support, α=-0.05, modifyA0=secant,
    nonlinearoptions=(maxiter=50, abstol=1e-5),
    copyOperators=false);
```

Both runs use `α = -0.05`, `Δt = 0.02` and `401` samples, including
the initial state, so the final time is `8.0`. The callbacks assemble fresh
matrices. With `copyOperators=false`, the incoming matrix is borrowed read-only;
these callbacks do not overwrite it or reuse previously returned matrix storage.

## 2. Newton-type tangent iteration

The analytical Total Lagrangian tangent contains a material and a geometric part:

```math
K_t=\int_{\Omega_0}G^T\bigl(B^TDB+G_\sigma\bigr)G\,b\,d\Omega_0.
```

To keep the *same internal force* as in the secant calculation, the right-hand-side callback applies

```math
b_{\mathrm{eq}}(u)=f_{\mathrm{ext}}+K_t(u)u-K_s(u)u.
```

**Discretization note:** the stress used in $G_\sigma$ is constructed as a LowLevelFEM field and interpolated to quadrature points. Thus this is an **analytical tangent / Newton-type iteration**, not a verified exact Jacobian of the discrete secant residual. An exactly consistent Gauss-point tangent requires evaluating stress from the deformation gradient at the same quadrature points (as in a dedicated element formulation).

```julia
function tangent(A, state)
    H, B = deformation(state)
    h11, h12, h21, h22 = H[1,1], H[1,2], H[2,1], H[2,2]
    e11 = h11 + (h11*h11 + h21*h21)/2
    e22 = h22 + (h12*h12 + h22*h22)/2
    γ12 = h12 + h21 + h11*h12 + h21*h22
    s11 = D2[1,1]*e11 + D2[1,2]*e22
    s22 = D2[2,1]*e11 + D2[2,2]*e22
    s12 = D2[3,3]*γ12
    Gσ = [s11 s12 0 0; s12 s22 0 0;
          0 0 s11 s12; 0 0 s12 s22]
    b*∫(Grad(U)⋅B'⋅D2⋅B⋅Grad(U) +
        Grad(U)⋅Gσ⋅Grad(U); gauss=4)
end

function equivalent_rhs(rhs, state)
    u = vec(state.fields[1].a)
    rhs.a[:,1] .+= state.operators.A0.A*u - secant(nothing, state).A*u
    return rhs
end
```

```julia
@time u_new, v_new = HHT(K, M, 0*u0, u0, v0, 401, 0.02;
    support, α=-0.05, modifyA0=tangent, modifyRHS=equivalent_rhs,
    nonlinearoptions=(maxiter=50, abstol=1e-5),
    copyOperators=false);
```

The matrix callback runs before the RHS callback, so `state.operators.A0`
contains the current tangent. The RHS argument is a private sample of the
prescribed load; the correction is applied to that sample on each evaluation.
Thus `K_t(u)u - b_eq(u) = K_s(u)u - f_ext`: both iterations target the same
internal-force residual. See the [callback contract](@ref transient-callbacks)
for the state members, update order and convergence controls.

## 3. Compare the responses

The two calculations use the same mass, HHT parameters, boundary conditions, initial conditions and internal-force residual. Compare the displacement histories; a smaller nonlinear iteration count alone does not establish greater accuracy.

```julia
println("Maximum DOF difference: ", maximum(abs.(u_sec.a - u_new.a)))
showDoFResults(u_sec; name="Cantilever — secant HHT", factor=1.0, visible=true)
showDoFResults(u_new; name="Cantilever — Newton-type HHT", factor=1.0, visible=true)
openPostProcessor()
```

The saved notebook run reports a maximum absolute displacement-DOF difference
of approximately `3.4e-7`. This is a comparison of the two iterations at the
chosen mesh, time step and tolerances. It does not establish time-step or mesh
convergence, or a generally valid performance advantage. The displayed `@time`
measurements include compilation and allocation effects.

The final `openPostProcessor()` call opens the Gmsh GUI. Omit that call when
running without a graphical display. Refine the mesh, halve the time increment
for the same final time, and tighten nonlinear tolerances to assess the response.
The interpolated-stress tangent qualification above remains essential when
interpreting Newton-type convergence.

## Related

- [Transient Linear Elasticity (2D)](transient-elasticity-2d.md)
- [Time Integration](../reference/transient.md)
- [Nonlinear Formulations with the DSL](nonlinear-large-deformation.md)
- [Linear Buckling and Prestressed Natural Frequencies](buckling-and-prestressed-modes.md)
