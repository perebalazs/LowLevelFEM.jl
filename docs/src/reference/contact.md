# Contact API

Contact geometry, weak-form operators and post-processing for slave-master contact problems.

## Overview

LowLevelFEM separates contact geometry and search from the numerical contact
formulation.

A `Contact` object describes one slave-master contact pair and stores the
geometry and search state required by contact integration. The actual weak-form
operator is created by `ContactGap(C)` and can be inserted directly into the
ordinary LowLevelFEM `∫` syntax.

The current configuration is

```math
x = X + u
```

on both the slave and master sides.

The contact formulation therefore follows the mathematical weak form directly.
For example, frictionless penalty contact is written as

```julia
G = ContactGap(C)

Kc = ∫(G ⋅ cn ⋅ G)
```

which corresponds to

```math
K_c
=
\int_{\Gamma_c}
G_q^T c_n G_q \, \mathrm d\Gamma .
```

The contact operator is evaluated directly at slave-side Gauss points. At every
integration point, LowLevelFEM computes the current slave position, performs a
closest-point projection onto the master manifold, evaluates the local contact
basis and constructs the slave-master kinematic operator.

No reduced `ContactVector` space is exposed in the public API.

Derived quantities such as penalty pressure, tangential traction, Lagrange
multipliers, stick/slip state or other constitutive contact results are not
stored as primary `Contact` data. They are constructed from the contact
kinematics and the chosen contact law.

---

# Constructing a contact pair

A contact pair is created with

```julia
C = contact(
    U;
    master="master",
    slave="slave",
    displacement=u
)
```

where:

- `U` is the displacement `Problem`,
- `master` is the master physical group,
- `slave` is the slave physical group,
- `u` is the current displacement field.

If `displacement` is omitted, a zero displacement field is used.

The convenience form

```julia
C = contact(
    u;
    master="master",
    slave="slave"
)
```

is also available.

The current geometry used by the contact search is always based on

```math
x = X + u.
```

---

## Main contact data

`Contact` is primarily a geometry and search object.

Useful public fields include:

| Field | Meaning |
| --- | --- |
| `master` | master physical group name |
| `slave` | slave physical group name |
| `U` | displacement `Problem` |
| `displacement` | displacement stored in the current contact state |
| `step` | current displacement step |
| `slave_nodes` | slave node tags |
| `master_element_tags` | nodal closest-point master element tags |
| `master_local_coordinates` | nodal master local coordinates |
| `master_points` | nodal projected master points |
| `gap` | nodal signed normal gap as a `ScalarField` |
| `gap_values` | nodal signed normal gap values |
| `n` | nodal contact normal as a `VectorField` |
| `t1` | first nodal tangent direction |
| `t2` | second nodal tangent direction in 3D |
| `active` | nodal active-contact mask |

The nodal fields above are refreshed by the full

```julia
updateContact!(C, u)
```

operation.

Gauss-point contact integration has its own projection data and warm-start cache.
It does not require the nodal gap fields to be recomputed at every nonlinear
iteration.

---

# Contact kinematics

At a slave integration point with local coordinate `ξ`, let

```math
x_s(\xi)
=
\sum_a N_a^s(\xi) x_a^s
```

be the current slave position.

Its closest master point has master local coordinate `η` and

```math
x_m(\eta)
=
\sum_b N_b^m(\eta) x_b^m .
```

The signed normal gap is

```math
g_n
=
(x_s-x_m)\cdot n.
```

With the default convention,

```math
g_n > 0
```

means open contact, while

```math
g_n < 0
```

means penetration.

For a frozen current contact geometry,

```math
\delta g_n
=
G_q \, \delta u_e .
```

In 3D, the full local contact basis may be written as

```math
Q_q
=
\begin{bmatrix}
n^T \\
t_1^T \\
t_2^T
\end{bmatrix},
```

and the full local contact operator is formed directly from the slave and
projected master interpolation:

```math
G_q
=
Q_q
\begin{bmatrix}
N_s I & -N_m I
\end{bmatrix}.
```

This operator is evaluated directly at Gauss points; no nodal contact matrix is
interpolated into the weak form.

---

# `ContactGap`

```julia
ContactGap(C; components=:normal, active=:current)
```

creates a contact kinematic operator for the LowLevelFEM weak-form DSL.

The default

```julia
G = ContactGap(C)
```

is a scalar normal-gap operator.

The full local relative-contact operator is requested with

```julia
G = ContactGap(C; components=:all)
```

with local component ordering

| Dimension | Ordering |
| ---: | --- |
| 2D | `(normal, tangent)` |
| 3D | `(normal, tangent1, tangent2)` |

The operator can be used directly in a matrix chain:

```julia
Kc = ∫(ContactGap(C) ⋅ cn ⋅ ContactGap(C))
```

or

```julia
G  = ContactGap(C; components=:all)
Dc = ContactStiffness(C, cn; ct=ct)

Kc = ∫(G ⋅ Dc ⋅ G)
```

The slave integration manifold is already stored in `C`, therefore `Γ="slave"`
is neither necessary nor accepted for `ContactGap` integration.

---

## Active Gauss points

The keyword

```julia
active=:current
```

is the default.

In this mode, a Gauss point contributes to the contact integral when its current
normal gap satisfies

```math
g_n \le \mathrm{activation\_tol}.
```

To integrate the entire slave candidate manifold, use

```julia
G = ContactGap(C; active=:all)
```

The Gauss-point active policy is evaluated during contact integration. It is
independent of the stored nodal mask `C.active`, which is only refreshed by a
full `updateContact!`.

---

# Penalty formulation

For frictionless penalty contact,

```math
\delta W_c
=
\int_{\Gamma_c}
\delta g_n \, c_n \, g_n \, \mathrm d\Gamma .
```

The corresponding tangent is

```math
K_c
=
\int_{\Gamma_c}
G_q^T c_n G_q \, \mathrm d\Gamma .
```

In LowLevelFEM:

```julia
G = ContactGap(C)

Kc = ∫(G ⋅ cn ⋅ G)
```

Let

```julia
r = nodePositionVector(U)
```

be the nodal reference-position field. Since `ContactGap` acts on the absolute
current position `r + u`, the contact residual is

```julia
rc = Kc * (r + u)
```

and the total residual is, for example,

```julia
R = K * u - f + rc
```

with tangent

```julia
A = K + Kc
```

for a frozen current contact geometry.

---

## Normal and tangential penalty stiffness

`ContactStiffness` creates a local contact constitutive coefficient without
requiring an explicit Julia matrix literal.

```julia
Dc = ContactStiffness(C, cn; ct=ct)
```

In 2D it represents

```math
\begin{bmatrix}
c_n & 0 \\
0   & c_t
\end{bmatrix},
```

and in 3D

```math
\begin{bmatrix}
c_n & 0   & 0 \\
0   & c_t & 0 \\
0   & 0   & c_t
\end{bmatrix}.
```

It is used with the full contact operator:

```julia
G = ContactGap(C; components=:all)
Kc = ∫(G ⋅ ContactStiffness(C, cn; ct=ct) ⋅ G)
```

Setting

```julia
ct = 0.0
```

removes tangential penalty stiffness.

A nonzero `ct` is a tangential penalty regularization. A physical Coulomb
friction law additionally requires tangential history and a stick/slip
algorithm.

---

# Updating the current configuration

There are two different update paths.

## Full update

```julia
updateContact!(C, u)
```

recomputes the complete nodal contact state, including:

- deformed contact geometry,
- nodal closest-point projections,
- nodal master element and local coordinates,
- nodal projected master points,
- nodal signed gaps,
- nodal normals and tangents,
- nodal active/inactive state.

The previous nodal master element and local coordinates are reused as a warm
start. The subsequent AABB search remains global and may select another master
element if a closer projection exists.

Use the full update when the nodal contact fields are needed, especially for
post-processing.

---

## Lightweight integration update with `updateFrom`

During nonlinear contact iteration, the Gauss-point weak form usually does not
need the complete nodal contact state.

The contact integral can therefore update only the deformed geometry required by
Gauss-point projection:

```julia
Kc = ∫(
    ContactGap(C) ⋅ cn ⋅ ContactGap(C);
    updateFrom=u
)
```

`updateFrom=u` performs a lightweight geometry update:

```text
u
↓
deformed nodal coordinates
↓
slave/master element coordinates
↓
master AABB tree
↓
Gauss-point projection and assembly
```

It deliberately does **not** recompute the nodal

```julia
C.gap
C.gap_values
C.n
C.t1
C.t2
C.active
```

fields.

Gauss-point closest-point projections are warm-started from the preceding
contact assembly.

Therefore, during a fast nonlinear iteration, do not interpret stored nodal
contact fields after using only `updateFrom=u`. Call

```julia
updateContact!(C, u)
```

when the full nodal state is required.

---

# Assembly options and performance

Contact bilinear forms use the same high-level assembly options as the ordinary
LowLevelFEM bilinear forms.

The default is direct CSC assembly:

```julia
Kc = ∫(
    G ⋅ cn ⋅ G;
    assembly=:csc,
    threads=:auto
)
```

The contact assembler uses:

- direct CSC storage,
- worker-local `nzval` buffers,
- parallel worker execution,
- parallel reduction,
- reusable Gauss-point closest-point warm starts,
- optional reusable CSC sparsity patterns.

The legacy triplet path remains available for validation:

```julia
Kc_ijv = ∫(
    G ⋅ cn ⋅ G;
    assembly=:ijv,
    threads=1
)
```

For debugging, the two paths can be compared numerically.

---

## CSC pattern reuse

The structural contact pattern can be reused between nonlinear iterations when
the current slave-master connectivity remains covered by the existing pattern.

First assemble normally:

```julia
Kc = ∫(
    G ⋅ cn ⋅ G;
    updateFrom=u,
    threads=:auto
)

pattern = copy(Kc.A)
```

Before every independent reuse, reset the numerical values:

```julia
fill!(pattern.nzval, 0.0)

Kc = ∫(
    G ⋅ cn ⋅ G;
    updateFrom=u,
    threads=:auto,
    csc_matrix=pattern
)
```

Assembly **adds** to the current `nzval` contents of a supplied pattern, so the
reset is required.

If the active slave-master connectivity changes in a way not represented by the
stored pattern, the contact assembler reports that the CSC pattern no longer
covers the current contact graph. Rebuild the pattern once by assembling without
`csc_matrix`, then reuse the new pattern.

The optional

```julia
element_chunk_size=:auto
```

keyword controls contact-element work partitioning for threaded assembly.

---

# Typical nonlinear penalty loop

A simple relaxed iteration can keep the nodal post-processing update out of the
inner loop:

```julia
G = ContactGap(C)

u_it = copy(u0)
ω = 0.5

Kc = ∫(
    G ⋅ cn ⋅ G;
    updateFrom=u_it,
    gauss=2,
    threads=:auto
)

pattern = copy(Kc.A)

for iter in 1:maxiter

    if iter > 1
        fill!(pattern.nzval, 0.0)

        Kc = ∫(
            G ⋅ cn ⋅ G;
            updateFrom=u_it,
            gauss=2,
            threads=:auto,
            csc_matrix=pattern
        )
    end

    rc = Kc * (r + u_it)
    R  = K * u_it - f + rc

    Δu = solveField(
        K + Kc,
        -R,
        support=support_increment
    )

    u_it = u_it + ω * Δu
end

u = u_it
```

After convergence, synchronize the complete nodal contact state once:

```julia
updateContact!(C, u)
```

This keeps the inner iteration focused on the Gauss-point contact weak form.

---

# Lagrange-multiplier coupling

The same normal-gap operator can be used in a mixed weak form.

For a scalar normal multiplier field `Λ`,

```julia
Gn = ContactGap(C)

B = ∫(
    Λ ⋅ Gn;
    updateFrom=u
)
```

corresponds to

```math
B
=
\int_{\Gamma_c}
N_\lambda^T G_n \, \mathrm d\Gamma .
```

The result is a rectangular `SystemMatrix` coupling the displacement trial space
to the multiplier test space.

The mixed contact assembler uses the same CSC, threading, Gauss-point projection
and warm-start infrastructure as penalty contact.

The active-set, complementarity or augmented-Lagrangian algorithm remains
separate from the geometric contact operator.

The legacy `LagrangeMultiplierField` keyword of `contact(...)` is retained only
for source compatibility; multiplier fields now enter the weak form directly.

---

# Several contact pairs

Several contact pairs can be stored in a `ContactSet`:

```julia
contacts = ContactSet(C1, C2, C3)
```

For penalty contact, each pair contributes its own weak-form matrix:

```julia
Kc =
    ∫(ContactGap(C1) ⋅ cn1 ⋅ ContactGap(C1); updateFrom=u) +
    ∫(ContactGap(C2) ⋅ cn2 ⋅ ContactGap(C2); updateFrom=u)
```

No reduced contact spaces need to be merged.

The full nodal state of every pair can be updated with

```julia
updateContact!(contacts, u)
```

when post-processing data are required.

---

# ContactGap post-processing

The same `ContactGap` object can be evaluated on a displacement field.

## Nodal evaluation

```julia
gap = ContactGap(C, u)
```

or equivalently

```julia
G = ContactGap(C)
gap = G(u)
```

returns a nodal `ScalarField` for the normal gap.

For the full local contact motion:

```julia
d = ContactGap(C, u; components=:all)
```

returns a nodal `VectorField` with local contact components.

The nodal evaluation uses the frozen nodal closest-point projections stored in
`C`. Therefore call

```julia
updateContact!(C, u)
```

first when the contact geometry should correspond to `u`.

---

## Gauss-point L2 projection

For post-processing based on the same slave Gauss-point contact kinematics used
by the weak form, specify `gauss`:

```julia
gap = ContactGap(C, u; gauss=2)
```

The gap is evaluated at slave Gauss points and globally L2-projected onto the
continuous slave-side Lagrange space.

The projected nodal coefficients satisfy

```math
\left(
\int_{\Gamma_s} N^T N \, \mathrm d\Gamma
\right)
g_h
=
\int_{\Gamma_s} N^T g_q \, \mathrm d\Gamma .
```

Available Gauss specifications follow the ordinary LowLevelFEM convention:

```julia
gauss = :full
gauss = :reduced
gauss = 0
gauss = 2
gauss = 8
```

Increasing the quadrature order improves the numerical projection of the
generally non-polynomial closest-point gap; it does not change the interpolation
order of the projected field.

For all local components:

```julia
d = ContactGap(
    C,
    u;
    components=:all,
    gauss=2
)
```

returns the globally projected local contact vector as a nodal `VectorField`.

---

## Penalty pressure

Penalty pressure does not require a separate contact-specific result type.

For the sign convention `g_n < 0` in penetration,

```julia
gap = ContactGap(C, u; gauss=2)

pressure = mapScalarField(
    g -> max(-cn * g, 0.0),
    gap
)
```

constructs the normal penalty pressure.

For tangential contact, the local tangential components returned by

```julia
ContactGap(C, u; components=:all, gauss=...)
```

can later be combined with the selected tangential constitutive or friction law
to obtain shear traction.

---

# Closest-point search

The master-side closest-point search uses an AABB tree.

The main search options are:

| Keyword | Default | Meaning |
| --- | ---: | --- |
| `aabb_padding` | `0.05` | relative expansion of master-element AABBs |
| `leaf_size` | `2` | maximum number of elements in an AABB leaf |
| `projection_tol` | `1e-10` | closest-point iteration tolerance |
| `projection_maxiter` | `40` | maximum projected Gauss-Newton iterations |

Basis information is cached for each master element type. Standard Lagrange
basis functions are represented locally by cached polynomial evaluators, which
avoids repeated Gmsh basis-function calls inside the closest-point iteration.

Gauss-point projection results also keep a warm-start state containing the
previous master element and master local coordinate.

---

# Topological stabilization near master boundaries

The ordinary closest-point path is preserved in element interiors.

Near a shared master vertex or edge, neighboring master elements may represent
essentially the same physical closest point. Small deformation changes can then
make the local representation switch repeatedly.

The contact search can stabilize projections near shared topological features:

```text
2D:
    element interior -> point-to-segment
    shared vertex    -> point-to-node

3D:
    face interior    -> point-to-face
    shared edge      -> point-to-edge
    shared vertex    -> point-to-node
```

The relevant options are:

| Keyword | Default | Meaning |
| --- | ---: | --- |
| `topology_tol` | `1e-3` | reference-space distance used to detect a shared edge or vertex |
| `topology_angle` | `45.0` | maximum incident-normal angle for treating the feature as smooth |

For example:

```julia
C = contact(
    U;
    master="master",
    slave="slave",
    displacement=u,
    topology_tol=1e-2
)
```

Setting

```julia
topology_tol=0.0
```

disables this stabilization.

The `topology_angle` criterion prevents smooth-feature stabilization from being
applied blindly across sharp geometric corners.

---

# Contact orientation

The signed gap uses the master-side normal:

```math
g_n = (x_s-x_m)\cdot n.
```

The orientation can be reversed with

```julia
normal_sign=-1.0
```

when required by the orientation of the master physical group.

For example:

```julia
C = contact(
    U;
    master="master",
    slave="slave",
    displacement=u,
    normal_sign=-1.0
)
```

---

# Self-contact

If

```julia
slave == master
```

self-contact mode is enabled automatically.

Local master elements connected to the slave point are excluded from the
closest-point search to prevent projection onto the point's own immediate
topological neighborhood.

The corresponding options are:

```julia
self_contact=true
self_exclusion_layers=1
```

Additional node-connected master-element layers can be excluded by increasing
`self_exclusion_layers`.

---

# Design principle

The contact API separates four layers:

```text
geometry and search
        ↓
Contact
        ↓
weak-form contact kinematics
        ↓
ContactGap(C)
        ↓
contact formulation
        ↓
penalty / multiplier / augmented / friction law
        ↓
post-processing
```

This keeps closest-point geometry independent of the numerical contact law.

The same `ContactGap` operator can therefore be reused in penalty,
Lagrange-multiplier, augmented and frictional formulations without changing the
underlying contact search.

---

# API reference

```@docs
Contact
ContactSet
ContactGap
ContactStiffness
contact
updateContact!
```
