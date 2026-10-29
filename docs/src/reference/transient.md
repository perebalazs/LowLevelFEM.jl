# Time Integration

```@meta
CurrentModule = LowLevelFEM
```

Use this page to choose a method, prepare data, find keyword defaults and
select a time increment. The complete function references remain at
[`HHT`](@ref), [`CDM`](@ref) and [`FDM`](@ref).

## Choose a method

| Equation | Method | Typical use |
|:--|:--|:--|
| `A2*ü + A1*u̇ + A0*u = rhs` | [`HHT`](@ref) | Implicit second-order dynamics; adjustable algorithmic damping. |
| `A2*ü + A1*u̇ + A0*u = rhs` | [`CDM`](@ref) | Explicit displacement updates; a stability-limited increment. |
| `A1*ẋ + A0*x = rhs` | [`FDM`](@ref) | First-order systems, including diffusion and heat conduction. |

The names K, M and C are positional aliases, not restrictions to elasticity:

| Generic operator | HHT / CDM argument | FDM argument |
|:--|:--|:--|
| A2: second derivative | M | Not present. |
| A1: first derivative | C, optional damping | C, required capacity. |
| A0: primary field | K | K |
| RHS | f | q |

```julia
u, v = HHT(K, M, C, f, u0, v0, n, Δt; support=bc, α=-0.1)
u, v = CDM(K, M, C, f, u0, v0, n, Δt;
           support=bc, dampingScheme=:central, velocityScheme=:central)
x = FDM(K, C, q, x0, n, Δt; support=bc, ϑ=0.5)
```

Omit positional C in HHT/CDM for zero physical damping, or provide it as
`damping=C` in that shorter calling form. The old positional BC-list forms
remain available; the keyword `support=bc` makes new calls easier to read.

## Inputs, time grid and results

| Input | Required layout |
|:--|:--|
| Operators | Compatible `SystemMatrix`s, optionally wrapped in `Symmetric`. |
| Loads | One constant column or exactly `n` columns in time order. |
| `U0`, `X0`, `V0` | Exactly one initial column; values initialize free coordinates. |
| `n` | Number of stored samples, including the initial one; at least 2. |
| `Δt`, `t0` | Constant positive increment and initial time (default `0.0`). |

The final time is `t0 + (n-1)*Δt`. Input Field time labels are **not**
interpolated; columns must already match this grid. After changing Δt, rebuild
sampled loads/supports and n to preserve the desired duration.

Loads and initial values accept nodal/element Fields, block `SystemVector`s
or field tuples in matrix order. Ordinary data use local solver components;
`Global(field)` or `Global(SystemVector(...))` marks Cartesian input. A local
load and `Global(...)` load can be added. Output consists of global nodal
Fields with one column per sample.

| Method | Single-field output | Multifield output |
|:--|:--|:--|
| HHT / CDM | `(u,v)` | `(u1,u2), (v1,v2)` in block order. |
| FDM | `x` | `(x1,x2)` in block order. |

## [Prescribed histories and constraints](@id transient-support)

### Supply each boundary condition once

`support` selects both the prescribed component and its history. Do not also
write that history into U0 or call `applyBoundaryConditions!(U0,bc)`.

```julia
bc = [BoundaryCondition("base"; field=U,
                        ux=(x,y,z,t)->0.001*sin(20t), uy=0.0)]
u, v = HHT(K, M, C, f, u0, v0, n, Δt; support=bc)
```

Accepted component values are:

| Value | Interpretation |
|:--|:--|
| Number | Constant in space and time. |
| `(x,y,z)->value` | Spatial function, constant in time. |
| `(x,y,z,t)->value` | Evaluated at each physical time sample. |
| `ScalarField` with 1 column | Constant component sampled over the mesh. |
| `ScalarField` with n columns | Complete sampled component history. |

If a function accepts both three and four arguments (optional arguments or
varargs), the spatial interpretation takes precedence. Use an explicit
four-argument closure for time dependence.

```julia
ux_history = ScalarField(U, "base", (x,y,z,t)->0.001*sin(20t);
                         steps=n, tmin=t0, tmax=t0+(n-1)*Δt)
bc = [BoundaryCondition("base"; field=U, ux=ux_history, uy=0.0)]
```

### Initial values and prescribed derivatives

- Boundary values override U0/X0 on constrained components.
- V0 initializes free first derivatives only. Prescribed velocities and
  accelerations follow the primary-variable history by finite differences.
- Interior stencils are centered; endpoint stencils are one-sided. Two samples
  define linear motion with zero acceleration; three use one common second
  difference; four or more permit second-order acceleration at the endpoints.
- To prescribe velocity data, integrate its scalar components with
  `integrate(v_history)` / `∫(v_history)` and use the results in a BC.

HHT integrates only the free coordinates with Newmark kinematics. Prescribed
inertia, damping and stiffness are first transferred to the effective load,
then receive the HHT endpoint weights. FDM includes the full capacity
increment `-C_fc*(d[j+1]-d[j])/Δt`.

### Coordinates, MPCs and approximation spaces

| Keyword | Default | Role |
|:--|:--|:--|
| `support` | `BoundaryCondition[]` | Prescribed component values/histories. |
| `mpc` | `MPC[]` | Kinematic ties between field DOFs. |
| `coordSys` | `CoordinateSystem[]` | Local nodal component bases. |
| `t0` | `0.0` | Initial physical time. |

The common transformation is `u_global=Q*T*u_reduced`: Q supplies coordinate
bases and T contains approximation/MPC mappings. MPCs identify local components;
prescribed slave values transfer to their masters at every sample. Each block
BC and MPC must identify `field=Problem`.

- Reduced-order and SEM interpolation use the existing cached transformations.
- Coordinates plus reduced/SEM interpolation on the same field are unsupported.
- SEM plus MPC on the same field is currently unsupported.
- Separate block fields may use these features independently.
- HHT/CDM need nonsingular mass on the remaining free coordinates; arbitrary
  mass-singular differential-algebraic systems are outside this interface.

Physical groups and components are resolved once per BC. Boundary-history
storage scales with prescribed rows times n, including transfer/restriction;
returned full solution histories still need nodal storage times n.

## HHT parameters

| Keyword | Default | Effect |
|:--|:--|:--|
| `α` | `0.0` | Algorithmic damping in `[-1/3,0]`; physical damping is in C. |
| `δ` | `nothing` | Optional override selecting `γ=0.5+δ`. |
| `γ` | `0.5-α` | Velocity parameter; a supplied value overrides the relation. |
| `β` | `(0.5+γ)^2/4` | Displacement parameter; can also be supplied explicitly. |
| `damping` | `nothing` | C by keyword, only when positional C is omitted. |

The standard family has `γ=0.5-α`, `β=(1-α)^2/4`. With `α=0` this is
average-acceleration Newmark. Explicit `δ=0` selects the historical independent
family. The implementation accepts custom finite `γ ≥ 0`, `β > 0`; those
choices require their own stability assessment.

## CDM scheme options

| Keyword | Default | Choices |
|:--|:--|:--|
| `dampingScheme` | `:lagged` | `:lagged` or `:central`. |
| `velocityScheme` | `:backward` | `:backward` or `:central`. |
| `damping` | `nothing` | C by keyword, only when positional C is omitted. |

| Choice | Meaning |
|:--|:--|
| Lagged damping | `C*(u[j]-u[j-1])/Δt`; first-order damping term, factor M. |
| Centered damping | `C*(u[j+1]-u[j-1])/(2Δt)`; second-order damping term, factor `M+Δt*C/2`. |
| Backward free velocity | Historical half-step estimate stored under the sample's time label. |
| Central free velocity | Centered sample-time derivative; one-sided final stencil; initial V0 retained. |

Constrained derivatives always follow the prescribed history, regardless of
`velocityScheme`. A diagonal mass uses division; otherwise CDM caches a
factorization, without automatic mass lumping. CDM has no linear-solver or
preconditioner keywords.

## FDM theta parameter

| Keyword | Default | Choices |
|:--|:--|:--|
| `ϑ` | `0.5` | Any finite value in `[0,1]`. |

| Value | Method | Accuracy and decay |
|:--|:--|:--|
| `0.0` | Explicit Euler | First order; stability-limited step. |
| `0.5` | Crank–Nicolson | Second order; weak damping of very stiff modes. |
| `1.0` | Backward Euler | First order; strong stiff-mode decay. |

For constant coefficients, FDM prepares `C_ff+ϑ*Δt*K_ff` once. With callbacks,
capacity is interpolated as `(1-ϑ)*A1_old+ϑ*A1_new`; each endpoint uses its own
`A0*x-rhs` contribution.

## [Choosing Δt: stability and accuracy](@id transient-timestep)

### CDM: respect the highest retained frequency

For a linear undamped positive-mass system, obtain the largest eigenvalue of
the **same constrained/reduced** pencil `K*ϕ=ω²*M*ϕ`. Then

```math
\Delta t < \frac{2}{\omega_{\max}}
          = \frac{T_{\min}}{\pi},
\qquad T_{\min}=\frac{2\pi}{\omega_{\max}}.
```

Lagged damping further limits the step. With classical modal damping and
`ξmax` bounding modal ratios over the retained spectrum, a conservative bound is

```math
\Delta t_{\mathrm{crit}}
= \frac{2}{\omega_{\max}}
  \left(\sqrt{1+\xi_{\max}^2}-\xi_{\max}\right).
```

Centered damping retains the undamped modal limit for nonnegative classical
damping. Use a safety factor below the relevant bound (for example 0.8) and
also check accuracy. These formulas are not a general stability certificate
for arbitrary nonsymmetric/nonproportionally damped operators.

```julia
# Undamped or centered classically damped CDM; use the run's own constraints.
λmax = largestEigenValue(K, M; support=bc, mpc=ties, coordSys=cs)
ωmax = sqrt(λmax)
Δt_limit = 2/ωmax
Δt = min(0.8*Δt_limit, T_interest/40)
```

Here `T_interest` is the shortest physical period that must be resolved,
including any faster imposed load or support motion. `largestEigenValue` uses
its optional ARPACK backend unless a direct solver is requested; load
`Arpack` first when needed. For lagged damping include its damping-dependent
factor in `Δt_limit`.

### HHT: unconditional linear stability still needs resolution

The standard HHT family is unconditionally stable for the usual linear system
with positive-definite mass, nonnegative stiffness and nonnegative physical
damping. This does not ensure accurate periods, peaks or wave propagation.

1. Start with roughly **20–50 samples per shortest period of interest**.
2. Resolve rapid load changes, support motion and expected contact events.
3. Inspect `:PE`, `:ADR` and `:SR` using [`HHTaccuracyAnalysis`](@ref).
4. Halve Δt, resample histories and compare the actual output quantities.

Algorithmic damping can suppress unwanted high-frequency modes; it cannot
recover a physical response that was not resolved. Nonlinear iteration also
has no unconditional convergence guarantee from the linear stability result.

### FDM: distinguish diffusion stability from decay accuracy

For `K*ϕ=λ*C*ϕ` with real nonnegative spectrum:

| Theta weight | Linear stability condition |
|:--|:--|
| `ϑ < 0.5` | `Δt < 2/((1-2ϑ)*λmax)`; leave a safety margin. |
| `ϑ ≥ 0.5` | No linear stability restriction on Δt for this spectral class. |

For complex spectra, examine the full theta stability region instead.
Even with Crank–Nicolson or backward Euler, resolve the shortest decay time
`1/λ` that matters and any faster forcing. Large `λ*Δt` can produce alternating
stiff modes with Crank–Nicolson. Use [`FDMaccuracyAnalysis`](@ref) and refine Δt
to check decay rates and peaks.

### When operators change

Reassess the CDM/explicit-FDM bound when callbacks change coefficients: contact
can increase the highest frequency abruptly. These APIs use **fixed Δt**, not
an automatic adaptive-step controller. A step that converges nonlinearly can
still be temporally inaccurate.

## [Linear solver and preconditioner options](@id transient-linear-solvers)

These four keywords apply to **HHT and FDM**, not CDM.

| Keyword | Default | Meaning |
|:--|:--|:--|
| `solver` | `:auto` | Backend in the table below. |
| `preconditioner` | `:auto` | `:auto`, `:none`, `:ich`, `:ilu`, or a preconstructed object. |
| `solveroptions` | `(;)` | Native iterative keywords, e.g. `(reltol=1e-9,maxiter=1000)`. |
| `preconditioneroptions` | `(;)` | Options used to construct `:ich` or `:ilu`. |

| Solver | Use |
|:--|:--|
| `:auto` | General sparse factorization; with symmetric declarations, Cholesky with LU fallback. |
| `:backslash` | Julia's generic factorization; diagonal systems use division. |
| `:lu` | LU factorization. |
| `:lu_no_ordering` | LU without the default column reordering. |
| `:cholesky` | Positive-definite symmetric system; no fallback. |
| `:ldlt` | Symmetric system; must be declared symmetric. |
| `:qr` | QR factorization. |
| `:cg` | Iterative conjugate gradients for a positive-definite symmetric system. |
| `:gmres` | Iterative GMRES for a general system. |

Wrap all contributing matrices in `Symmetric(...)` for a symmetric effective
system: K/M and supplied C in HHT, K/C in FDM. Upper- or lower-triangle storage
is supported. A symmetric declaration does not prove positive definiteness.

| Preconditioner | Dependency / behavior |
|:--|:--|
| `:auto` | CG selects `:ich`; GMRES selects `:ilu`. |
| `:none` | No preconditioner. |
| `:ich` | `Preconditioners`; constructor `memory` defaults to 2. |
| `:ilu` | `IncompleteLU`; forwards constructor options such as `τ`. |
| Supplied object | Reused directly; constructor options must be empty. |

Both iterative solvers need `using IterativeSolvers`. Load the selected
preconditioner package explicitly. Automatic CG preconditioning requires the
symmetric declaration. Nonempty option tuples and actual preconditioners are
rejected with direct solvers. `preconditioneroptions` also cannot accompany
`:none` or a preconstructed object. Native `log`, `Pl`, `Pr`, `initially_zero`
are managed internally; other supported native options depend on the backend.

```julia
using IterativeSolvers, Preconditioners
u, v = HHT(Symmetric(K), Symmetric(M), Symmetric(C),
           f, u0, v0, n, Δt;
           solver=:cg, solveroptions=(reltol=1e-10, maxiter=2000),
           preconditioner=:ich, preconditioneroptions=(memory=2,))
```

Fixed operators reuse factors/preconditioners. HHT also prepares a separate
initial mass solve with the same settings unless the mass is diagonal.
See [`solveField`](@ref) for the shared backend contract.

## [Operator updates and nonlinear iteration](@id transient-callbacks)

### Callback options

All callbacks default to `nothing` and use `callback(operator,state)`.
Extra constants, contact objects and model parameters can be captured by a
closure or callable object; no `Function` type restriction is imposed.

| Keyword | HHT / CDM | FDM |
|:--|:--|:--|
| `modifyA2` | Update M, the second-derivative operator. | Only `nothing` is accepted. |
| `modifyA1` | Update C, the first-derivative operator. | Update C, the capacity operator. |
| `modifyA0` | Update K, the primary-field operator/tangent. | Update K, the primary-field operator/tangent. |
| `modifyRHS` | Update f at the evaluated state. | Update q at the trial endpoint. |
| `nonlinearoptions` | Residual iteration settings, default `(;)`. | Same settings. |

With no callbacks, keep `nonlinearoptions=(;)` for the fast linear path.
HHT/FDM iterate within the time step. CDM retains explicit displacement
updates and only iterates acceleration/centered-velocity consistency when
needed. Changing coefficients does not remove explicit stability limits.

### What is in `state`?

| Member | Meaning |
|:--|:--|
| `fields` | Tuple of trial global nodal Fields, even for one field. |
| `velocity` | Tuple of global first derivatives; `nothing` at FDM initialization. |
| `acceleration` | Tuple of global second derivatives for HHT/CDM; `nothing` for FDM. |
| `problems` | Tuple of corresponding Problem objects, including materials. |
| `t` | Current physical evaluation time. |
| `dt` | Constant time increment. |
| `step` | Zero-based sample index; the initial sample is 0. |
| `iteration` | 0 for the predictor, then the completed correction count. |
| `phase` | `:initial`, `:step` for HHT/FDM, or `:explicit` for CDM. |
| `previous` | Previous accepted field snapshot; `nothing` at initialization. |
| `operators` | Current full/global A2, A0, A1; exposed as updates proceed. |

Callbacks run in order **A2, A0, A1, RHS**. For example, the A1 callback can use
updated A2/A0 to build Rayleigh damping.

### Inputs, returns and accepted state

| Callback | Input | Return |
|:--|:--|:--|
| Matrix callback | Private copy of the last accepted full operator. | `SystemMatrix`, optionally `Symmetric`, or `nothing`. |
| RHS callback | Private one-column original prescribed RHS sample at `state.t`. | Global one-column Field, tuple, `SystemVector`, or `nothing`. |

Return the private matrix argument to retain in-place edits. Returning
`nothing` keeps the accepted matrix or the original prescribed RHS sample.
Every trial starts from the accepted base; changes do not accumulate across
rejected/unequilibrated trials.

Callback data use **global Cartesian components**, including RHS values.
Do not preapply coordinates, MPCs, Dirichlet elimination or basis reduction.
Treat state Fields, Problems and `state.operators` as read-only. Callbacks may
run repeatedly at the same time; do not commit irreversible material-history
updates on the assumption that each call is an accepted step.

Numerical values and sparse connectivity may change. Problem identities,
field order, dimensions, mesh connectivity, supports, MPCs, coordinates and
approximation spaces stay fixed. Previously removed auxiliary MPC DOFs cannot
be activated. Equal operators reuse their cached projection/factorization;
RHS-only updates do not invalidate a matrix factorization.

### Nonlinear controls

| `nonlinearoptions` member | Default | Meaning |
|:--|:--|:--|
| `maxiter` | `30` | Maximum corrections; positive integer. |
| `reltol` | `1e-9` | Nonnegative relative residual tolerance. |
| `abstol` | `1e-11` | Nonnegative absolute tolerance in equilibrium-equation units. |
| `relaxation` | `1.0` | Correction fraction in `(0,1]`. |

The acceptance condition is `norm(residual) ≤ abstol + reltol*scale`, using
assembled equilibrium terms for the scale. A small correction alone does not
establish equilibrium. A failure throws the method, step, time and residual;
the solution is not advanced on an unconverged trial.

### Newton/Picard example

For internal force `g(u)`, return its tangent J(u) as A0 and
`f_external-g(u)+J(u)*u` as RHS. Picard coefficient updates are also possible.
State-dependent A2/A1 require their own consistent linearization for Newton
convergence; the residual loop does not construct those derivatives.

```julia
# Scalar cubic internal force: k*u + λ*u^3; k, λ and P are captured.
updateA0 = (A, state) -> begin
    u = vec(state.fields[1].a)
    SystemMatrix(spdiagm(0 => k .+ 3λ .* u.^2), P)
end
updateRHS = (rhs, state) -> begin
    out = copy(rhs)
    out.a .+= 2λ .* state.fields[1].a.^3
    out
end
u, v = HHT(K, M, C, f, u0, v0, n, Δt;
           support=bc, modifyA0=updateA0, modifyRHS=updateRHS,
           nonlinearoptions=(reltol=1e-10, maxiter=30))
```

## Multifield remote translation and rotation

Loads and initial data follow block order, for example
`(u1,u2), (v1,v2) = HHT(Kblock,Mblock,Cblock,F,U0,V0,n,Δt; ...)`.
FDM returns `(x1,x2)` and uses the same preparation.

For physical displacement `u_phys=u+R*φ`, apply the same rigid map to **every**
operator and load:

```julia
R = rigidRotationMap(mpc_u, mpc_φ)
block(A) = SystemMatrix([A A*R; R'*A R'*A*R])
Kblock, Mblock, Cblock = block(Ku), block(Mu), block(Cu)
F = SystemVector([fu, R'*fu])

(u,φ), (v,ω) = HHT(Kblock, Mblock, Cblock, F, U0, V0, n, Δt;
                   support=bc, mpc=[mpc_u,mpc_φ])
u_phys = u + R*φ
v_phys = v + R*ω
```

Transforming stiffness alone omits rotational inertia/coupling. No artificial
independent mass is needed on the auxiliary rotation field: its inertia comes
from `R'*Mu*R`. Inactive auxiliary DOFs are removed before time integration.

## Accuracy-analysis helpers

| Function | Spectrum | Available analysis |
|:--|:--|:--|
| [`CDMaccuracyAnalysis`](@ref) | `K*ϕ=ω²*M*ϕ` | `:SR`, `:PDR`, `:ADR`, `:PE`; match `dampingScheme`. |
| [`HHTaccuracyAnalysis`](@ref) | `K*ϕ=ω²*M*ϕ` | Same indicators; includes modal and Rayleigh damping. |
| [`FDMaccuracyAnalysis`](@ref) | `K*ϕ=λ*C*ϕ` | `type=:SR` for the theta method. |

Each returns plotting vectors `(x,y)` and accepts `n=100` sample points.
For HHT/CDM `:SR`, x is the **natural logarithm** of `ω*Δt/(2π)`;
other HHT/CDM curves use `ω*Δt/(2π)`. FDM uses `λ*Δt`.

```julia
x, radius = HHTaccuracyAnalysis(ωmin, ωmax, Δt, :SR;
                               α=-0.1, ξ=0.02,
                               massDamping=0.0, stiffnessDamping=1e-4)
```

Use the same damping parameters as the integrator. These scalar analyses assume
modal/proportional damping; a general nonproportional matrix cannot be inferred
from the undamped frequency interval alone. The helper docstrings list every
parameter and the CDM helper's legacy mass-damping convention.

See also the [2D transient elasticity tutorial](../tutorials/transient-elasticity-2d.md)
and the shared [Solution API](solution.md).
