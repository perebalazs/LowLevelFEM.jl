# Contact API

Contact kinematics and algebraic interface for node-to-manifold contact problems.

## Overview

LowLevelFEM separates contact geometry and kinematics from the numerical contact
method. A `Contact` object describes one slave-master contact pair in the current
configuration and provides the algebraic operators required by penalty,
Lagrange-multiplier and augmented formulations.

The current geometry is evaluated as

```math
x = X + u
```

on both the slave and master sides.

The contact search determines, for every slave contact node:

- the closest master point,
- the signed normal gap,
- the local normal and tangential basis,
- the active/inactive state,
- the reduced contact-space kinematic operator.

No contact pressure, traction, multiplier solution, stick/slip state or other
derived contact result is stored in `Contact`. These quantities are constructed
explicitly from the supplied operators in the calling code.

---

# Reduced contact space

Let

```math
V_u
```

denote the global displacement space and

```math
V_c
```

the reduced local contact space.

For `nc` slave contact nodes and spatial dimension `pdim`, the reduced contact
space has

```math
n_c = nc \, pdim
```

components.

The local ordering is:

| Dimension | Per-node ordering |
| ---------:| ----------------- |
| 2D        | `[normal, tangent]` |
| 3D        | `[normal, tangent1, tangent2]` |

Thus the complete ordering is

```text
2D: [n1, t1, n2, t2, ...]
3D: [n1, t11, t21, n2, t12, t22, ...]
```

The normal component is always the first component belonging to a contact node.

---

## `ContactVector`

`ContactVector` is an algebraic vector living in the reduced contact space.

It is intentionally different from `ScalarField`, `VectorField` and
`TensorField`, because its size is determined by the active contact
discretization rather than by the nodal finite-element field layout.

For example:

```julia
g = contact.g
```

returns the reduced contact gap vector.

The basic contact-space operations include:

```julia
g1 + g2
g1 - g2
α * g
g / α
norm(g)
dot(g1, g2)
```

A `ContactVector` cannot be added directly to a finite-element field. It must
first be mapped through an appropriate `SystemMatrix`.

Typical mappings are:

```math
C : V_c \rightarrow V_c
```

and

```math
G^T : V_c \rightarrow V_u.
```

Therefore

```julia
C * g
```

returns another `ContactVector`, while

```julia
G' * C * g
```

returns a nodal displacement-space `VectorField`.

---

# Contact kinematics

For one contact pair, the kinematic operator is

```math
G : V_u \rightarrow V_c.
```

If `ndofs(U)` is the number of displacement degrees of freedom, then

```math
G \in \mathbb{R}^{(nc\,pdim)\times ndofs(U)}.
```

The operator maps a global displacement increment to local relative contact
motion.

In 2D:

```math
G \, \Delta u
=
\begin{bmatrix}
\Delta u_n \\
\Delta u_t
\end{bmatrix}
```

for each slave contact node.

In 3D:

```math
G \, \Delta u
=
\begin{bmatrix}
\Delta u_n \\
\Delta u_{t_1} \\
\Delta u_{t_2}
\end{bmatrix}.
```

The signed normal gap is

```math
g_n = (x_s - x_m)\cdot n.
```

With the default convention:

```math
g_n > 0
```

means an open contact point, while

```math
g_n < 0
```

means penetration.

---

# Constructing a contact pair

A contact pair is created with

```julia
contact_pair = contact(
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
contact_pair = contact(
    u;
    master="master",
    slave="slave"
)
```

is also available.

---

## Main contact data

The most important fields of a `Contact` object are:

| Field | Meaning |
| ----- | ------- |
| `gap` | signed normal gap as a `ScalarField` |
| `gap_values` | signed normal gap values indexed by contact node |
| `g` | reduced local gap as a `ContactVector` |
| `G` | kinematic operator `Vu -> Vc` |
| `C` | local contact-space operator `Vc -> Vc` |
| `E` | optional embedding `Vc -> Vλ` for multiplier contact |
| `n` | contact normal as a `VectorField` |
| `t1` | first tangent direction |
| `t2` | second tangent direction in 3D |
| `active` | active contact-point mask |
| `slave_nodes` | slave node tags |
| `master_element_tags` | projected master element tags |
| `master_local_coordinates` | master local coordinates of the projections |
| `master_points` | projected master points in physical space |

The contact object also stores the normal and tangential stiffness definitions
`cn` and `ct`, together with their evaluated nodal values.

---

# Contact-space operator `C`

The local operator `C` is constructed from the normal and tangential contact
stiffness values.

For frictionless penalty contact, use

```julia
cn = 1.0e6
ct = 0.0
```

or any appropriate problem-dependent stiffness.

Both `cn` and `ct` may be:

- a scalar number,
- a `ScalarField`,
- a function `f(x,y,z)`.

For isotropic tangential regularization, the same `ct` value is used in each
local tangential direction.

---

# Penalty formulation

For a penalty formulation, the contact contribution can be written directly
with the operators supplied by `Contact`.

Let

```math
g \in V_c,
```

```math
C : V_c \rightarrow V_c,
```

and

```math
G : V_u \rightarrow V_c.
```

The local contact traction-like vector is

```math
p = -C g.
```

The corresponding global contact residual contribution is

```math
r_c = -G^T p = G^T C g.
```

The contact tangent is

```math
K_c = G^T C G.
```

In LowLevelFEM:

```julia
(; G, C, g) = contact_pair

p  = -C * g
rc = -G' * p
Kc = G' * C * G

r = K * u - f + rc
A = K + Kc
```

The algebra therefore follows the mathematical formulation directly, without a
separate penalty-contact solver wrapper.

---

# Updating the current contact configuration

During a nonlinear iteration, update the contact geometry with

```julia
updateContact!(contact_pair, u)
```

This recomputes:

- closest-point projections,
- contact normals and tangents,
- signed gaps,
- active contact points,
- `g`,
- `G`,
- `C`,
- optional multiplier embedding `E`.

The previous master element and local coordinates are reused as a warm start.
The subsequent AABB search remains global and may select another master element
if a closer projection is found.

A typical nonlinear penalty loop therefore contains:

```julia
updateContact!(contact_pair, u_it)

(; G, C, g) = contact_pair

Kc = G' * C * G
rc = G' * C * g

r = K * u_it - f + rc
A = K + Kc
```

---

# Lagrange-multiplier contact

A vector-valued multiplier `Problem` can be associated with the contact pair:

```julia
contact_pair = contact(
    U;
    master="master",
    slave="slave",
    displacement=u,
    LagrangeMultiplierField=Λ
)
```

The multiplier field must have the same local component dimension as the
contact space:

- 2 components in 2D,
- 3 components in 3D.

The contact kinematics remain unchanged:

```math
G : V_u \rightarrow V_c.
```

An additional embedding operator is provided:

```math
E : V_c \rightarrow V_\lambda.
```

This maps the reduced contact quantities into the finite-element multiplier
field.

The multiplier coupling matrix is therefore

```math
B = E G,
```

and the multiplier-space gap residual is

```math
g_\lambda = E g.
```

In LowLevelFEM:

```julia
(; G, g, E) = contact_pair

B  = E * G
gλ = E * g
```

These operators can then be used directly in a multifield block system.

For example, the linearized mixed system has the structure

```math
\begin{bmatrix}
K & B^T \\
B & 0
\end{bmatrix}
\begin{bmatrix}
\Delta u \\
\Delta \lambda
\end{bmatrix}
=
-
\begin{bmatrix}
r_u \\
r_\lambda
\end{bmatrix}.
```

The precise active-set or complementarity algorithm is intentionally not part
of `Contact`.

---

# Several contact pairs

Several independent slave-master pairs can be stored in a `ContactSet`:

```julia
contacts = ContactSet(c1, c2, c3)
```

Each `Contact` retains its own reduced contact space and its own `G`, `C`, `g`
and optional multiplier embedding `E`.

The reduced spaces are not merged automatically.

This is particularly important for multifield multiplier formulations, where
different contact pairs may use different multiplier fields.

For penalty contact, the global contributions may be assembled directly:

```julia
Kc = sum(c.G' * c.C * c.G for c in contacts)
rc = sum(c.G' * c.C * c.g for c in contacts)
```

All pairs can be updated with the same displacement field:

```julia
updateContact!(contacts, u)
```

For multiplier contact, the coupling operators remain pair-specific:

```julia
B1 = contacts[1].E * contacts[1].G
B2 = contacts[2].E * contacts[2].G
```

This preserves the `model` and `test_model` metadata required by the multifield
assembly.

---

# Closest-point search

The master-side closest-point search uses an AABB tree.

The main search options are:

| Keyword | Default | Meaning |
| ------- | ------: | ------- |
| `aabb_padding` | `0.05` | relative expansion of master-element AABBs |
| `leaf_size` | `2` | maximum number of elements in an AABB leaf |
| `projection_tol` | `1e-10` | closest-point iteration tolerance |
| `projection_maxiter` | `40` | maximum projected Gauss-Newton iterations |

Basis information is cached for each master element type. The iterative
closest-point search uses the cached local polynomial evaluator and reusable
workspaces, avoiding repeated Gmsh basis-function calls inside the projection
loop.

---

# Topological stabilization near master boundaries

The ordinary node-to-manifold closest-point path is preserved in element
interiors.

Near a shared master vertex or edge, however, two neighboring master elements
may represent essentially the same physical closest point. Small changes in
the deformation can then make the projection switch repeatedly between the two
elements.

To reduce this ambiguity, the contact search can replace the ordinary local
representation near shared topological features by a stable representation:

```text
2D:
    element interior -> node-to-segment
    shared vertex    -> node-to-node

3D:
    face interior    -> node-to-face
    shared edge      -> node-to-edge
    shared vertex    -> node-to-node
```

The relevant geometry options are:

| Keyword | Default | Meaning |
| ------- | ------: | ------- |
| `topology_tol` | `1e-3` | reference-space distance used to detect proximity to a shared edge or vertex |
| `topology_angle` | `45.0` | maximum angle in degrees between incident normals for treating the feature as smooth |

For example:

```julia
contact_pair = contact(
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

disables this stabilization and leaves the ordinary closest-point path
unchanged.

The `topology_angle` criterion prevents smooth-feature stabilization from
being applied blindly across sharp geometric corners.

---

# Active contact points

A contact point is marked active when

```math
g_n \le \text{activation\_tol}.
```

The default is

```julia
activation_tol = 0.0
```

which corresponds to geometric penetration or exact contact.

A positive activation tolerance may be useful when an algorithm should include
points that are still separated by a small distance.

The active-state information is available as

```julia
contact_pair.active
```

while the signed normal gaps are available directly as

```julia
contact_pair.gap_values
```

---

# Contact orientation

The signed gap uses the master-side normal:

```math
g_n = (x_s - x_m)\cdot n.
```

The orientation can be reversed with

```julia
normal_sign = -1.0
```

when required by the orientation of the master physical group.

For example:

```julia
contact_pair = contact(
    U;
    master="master",
    slave="slave",
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
closest-point search to prevent the node from projecting onto its own immediate
topological neighborhood.

The corresponding options are:

```julia
self_contact = true
self_exclusion_layers = 1
```

Additional node-connected master-element layers can be excluded by increasing
`self_exclusion_layers`.

---

# Contact quantities and post-processing

`Contact` stores only primitive geometric and algebraic contact information.

Derived quantities should be constructed explicitly from the contact method.

For penalty contact, for example:

```julia
p = -contact_pair.C * contact_pair.g
```

For multiplier contact, the multiplier field itself is an independent unknown
of the mixed problem.

The contact normal and tangent fields are available as:

```julia
contact_pair.n
contact_pair.t1
contact_pair.t2
```

and the scalar normal gap field as:

```julia
contact_pair.gap
```

Ordinary LowLevelFEM post-processing functions can therefore be used after the
desired contact quantity has been mapped to a standard finite-element field.

---

# Design principle

The `Contact` interface deliberately separates three layers:

```text
geometry / kinematics
        ↓
G, g, n, t1, t2, active
        ↓
contact formulation
        ↓
penalty / Lagrange / augmented / friction law
        ↓
derived results
```

This keeps the geometric contact search independent of the numerical contact
method.

The same contact kinematics can therefore be reused by several formulations
without changing the closest-point algorithm or the finite-element model.

---

# API reference

```@docs
Contact
ContactSet
ContactVector
contact
updateContact!
```
