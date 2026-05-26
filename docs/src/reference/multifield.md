# Multifield API

Weak-form DSL and multifield assembly interface.

## Operator size and component ordering

This section summarizes the output size and component ordering of the basic
differential operators used by the weak-form DSL.

All operator vectors are written in the order used internally by LowLevelFEM.

---

## Scalar fields

Let

```math
p = p(x,y,z)
```

be a scalar field.

| Dimension | Operator  | Size | Component ordering |
| ---------:| --------- | ----:| ------------------ |
| 1D        | `Grad(P)` | `1`  | `[p,x]`            |
| 2D        | `Grad(P)` | `2`  | `[p,x, p,y]`       |
| 3D        | `Grad(P)` | `3`  | `[p,x, p,y, p,z]`  |
| 1D/2D/3D  | `Id(P)`   | `1`  | `[p]`              |

`Div(P)`, `Curl(P)` and `SymGrad(P)` are not defined for scalar fields.

---

## Vector fields

Let

```math
u =
\begin{bmatrix}
u_x \\
u_y \\
u_z
\end{bmatrix}
```

be a vector field. In 2D, only `ux` and `uy` are present.

---

## `Grad(Pu)` for vector fields

`Grad(Pu)` returns the full displacement gradient in component-major ordering.

### 2D vector field

```math
u =
\begin{bmatrix}
u_x \\
u_y
\end{bmatrix}
```

| Operator   | Size | Component ordering         |
| ---------- | ----:| -------------------------- |
| `Grad(Pu)` | `4`  | `[ux,x, ux,y, uy,x, uy,y]` |

Equivalent matrix form:

```math
\nabla u =
\begin{bmatrix}
u_{x,x} & u_{x,y} \\
u_{y,x} & u_{y,y}
\end{bmatrix}
```

Flattened as:

```math
[u_{x,x}, u_{x,y}, u_{y,x}, u_{y,y}]
```

### 3D vector field

| Operator   | Size | Component ordering                                       |
| ---------- | ----:| -------------------------------------------------------- |
| `Grad(Pu)` | `9`  | `[ux,x, ux,y, ux,z, uy,x, uy,y, uy,z, uz,x, uz,y, uz,z]` |

Equivalent matrix form:

```math
\nabla u =
\begin{bmatrix}
u_{x,x} & u_{x,y} & u_{x,z} \\
u_{y,x} & u_{y,y} & u_{y,z} \\
u_{z,x} & u_{z,y} & u_{z,z}
\end{bmatrix}
```

Flattened as:

```math
[u_{x,x}, u_{x,y}, u_{x,z},
 u_{y,x}, u_{y,y}, u_{y,z},
 u_{z,x}, u_{z,y}, u_{z,z}]
```

---

## `SymGrad(Pu)` for vector fields

`SymGrad(Pu)` returns the engineering strain vector.

### 2D vector field

| Operator      | Size | Component ordering |
| ------------- | ----:| ------------------ |
| `SymGrad(Pu)` | `3`  | `[εxx, εyy, γxy]`  |

Explicitly:

```math
\operatorname{SymGrad}(u) =
\begin{bmatrix}
u_{x,x} \\
u_{y,y} \\
u_{x,y} + u_{y,x}
\end{bmatrix}
```

### 3D vector field

| Operator      | Size | Component ordering               |
| ------------- | ----:| -------------------------------- |
| `SymGrad(Pu)` | `6`  | `[εxx, εyy, εzz, γxy, γyz, γzx]` |

Explicitly:

```math
\operatorname{SymGrad}(u) =
\begin{bmatrix}
u_{x,x} \\
u_{y,y} \\
u_{z,z} \\
u_{x,y} + u_{y,x} \\
u_{y,z} + u_{z,y} \\
u_{z,x} + u_{x,z}
\end{bmatrix}
```

The shear components are engineering shear strains.

---

## `Div(Pu)` for vector fields

| Dimension | Operator  | Size | Component ordering     |
| ---------:| --------- | ----:| ---------------------- |
| 2D        | `Div(Pu)` | `1`  | `[ux,x + uy,y]`        |
| 3D        | `Div(Pu)` | `1`  | `[ux,x + uy,y + uz,z]` |

---

## `Curl(Pu)` for vector fields

### 2D vector field

| Operator   | Size | Component ordering |
| ---------- | ----:| ------------------ |
| `Curl(Pu)` | `1`  | `[uy,x - ux,y]`    |

### 3D vector field

| Operator   | Size | Component ordering                        |
| ---------- | ----:| ----------------------------------------- |
| `Curl(Pu)` | `3`  | `[uz,y - uy,z, ux,z - uz,x, uy,x - ux,y]` |

That is:

```math
\nabla \times u =
\begin{bmatrix}
u_{z,y} - u_{y,z} \\
u_{x,z} - u_{z,x} \\
u_{y,x} - u_{x,y}
\end{bmatrix}
```

---

## Tensor fields

A second-order tensor field is stored as a full `3×3` tensor, even in many
2D workflows.

The internal tensor component ordering is column-major:

```math
T =
\begin{bmatrix}
T_{11} & T_{12} & T_{13} \\
T_{21} & T_{22} & T_{23} \\
T_{31} & T_{32} & T_{33}
\end{bmatrix}
```

stored as:

```math
[T_{11}, T_{21}, T_{31},
 T_{12}, T_{22}, T_{32},
 T_{13}, T_{23}, T_{33}]
```

| Stored index | Tensor component |
| ------------:| ---------------- |
| `1`          | `T11`            |
| `2`          | `T21`            |
| `3`          | `T31`            |
| `4`          | `T12`            |
| `5`          | `T22`            |
| `6`          | `T32`            |
| `7`          | `T13`            |
| `8`          | `T23`            |
| `9`          | `T33`            |

---

## `TensorDiv(P)` for tensor fields

Let `T` be a second-order tensor field.

| Dimension | Operator       | Size | Component ordering                                                      |
| ---------:| -------------- | ----:| ----------------------------------------------------------------------- |
| 2D        | `TensorDiv(P)` | `2`  | `[T11,x + T12,y, T21,x + T22,y]`                                        |
| 3D        | `TensorDiv(P)` | `3`  | `[T11,x + T12,y + T13,z, T21,x + T22,y + T23,z, T31,x + T32,y + T33,z]` |

In index notation:

```math
(\operatorname{div} T)_i =
\frac{\partial T_{ij}}{\partial x_j}
```

---

## Voigt convention

LowLevelFEM uses the following 3D Voigt ordering for symmetric second-order
tensors:

```math
[xx, yy, zz, xy, yz, zx]
```

That is:

| Voigt index | Component |
| -----------:| --------- |
| `1`         | `xx`      |
| `2`         | `yy`      |
| `3`         | `zz`      |
| `4`         | `xy`      |
| `5`         | `yz`      |
| `6`         | `zx`      |

For stress-like tensors:

```math
[S_{xx}, S_{yy}, S_{zz}, S_{xy}, S_{yz}, S_{zx}]
```

For engineering strain-like vectors:

```math
[\varepsilon_{xx}, \varepsilon_{yy}, \varepsilon_{zz},
 \gamma_{xy}, \gamma_{yz}, \gamma_{zx}]
```

where

```math
\gamma_{xy} = 2\varepsilon_{xy}
```

and similarly for the other shear components.

---

## Useful nonlinear 2D mapping

For large-displacement 2D formulations embedded in 3D tensor notation,
`Grad(Pu)` has four components:

```math
[u_{x,x}, u_{x,y}, u_{y,x}, u_{y,y}]
```

The variation of the Green-Lagrange strain can be written as

```math
\delta E_\mathrm{voigt} = A(F) \, \nabla \delta u
```

with the 3D Voigt ordering

```math
[xx, yy, zz, xy, yz, zx]
```

and

```math
F =
\begin{bmatrix}
F_{11} & F_{12} & 0 \\
F_{21} & F_{22} & 0 \\
0      & 0      & F_{33}
\end{bmatrix}.
```

Then

```math
A(F) =
\begin{bmatrix}
F_{11} & 0      & F_{21} & 0      \\
0      & F_{12} & 0      & F_{22} \\
0      & 0      & 0      & 0      \\
F_{12} & F_{11} & F_{22} & F_{21} \\
0      & 0      & 0      & 0      \\
0      & 0      & 0      & 0
\end{bmatrix}.
```

Thus the material tangent contribution can be assembled as:

```julia
Kmat = ∫(Grad(Pu) ⋅ A' ⋅ D ⋅ A ⋅ Grad(Pu); Ω="body")
```

with dimensions:

| Quantity   | Size  |
| ---------- | -----:|
| `Grad(Pu)` | `4`   |
| `A`        | `6×4` |
| `D`        | `6×6` |
| `A'`       | `4×6` |

---

## 2D restriction matrix

For extracting the in-plane components `[xx, yy, xy]` from the 3D Voigt vector
`[xx, yy, zz, xy, yz, zx]`, use:

```math
R =
\begin{bmatrix}
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 0 \\
0 & 0 & 1 \\
0 & 0 & 0 \\
0 & 0 & 0
\end{bmatrix}.
```

Then:

```math
D_{2D} = R^T D R
```

and

```math
S_{2D} = R^T S_\mathrm{voigt}.
```

The corresponding 2D material tangent contribution is:

```julia
A2 = R' * A
D2 = R' * D * R

Kmat = ∫(Grad(Pu) ⋅ A2' ⋅ D2 ⋅ A2 ⋅ Grad(Pu); Ω="body")
```

# Embedded surface operators

## `SurfaceGrad(P)` for scalar fields

Let

```math
p = p(x,y,z)
```

be a scalar field defined on a surface embedded in 3D.

`SurfaceGrad(P)` returns the tangential surface gradient.

### Embedded surface in 3D

| Operator         | Size | Component ordering |
| ---------------- | ----:| ------------------ |
| `SurfaceGrad(P)` | `2`  | `[p,1, p,2]`       |

where:

* `1` and `2` denote the local tangent directions
  `(t₁,t₂)` evaluated at the Gauss point.

Mathematically:

```math
\nabla_\Gamma p
=
\begin{bmatrix}
\partial p/\partial s_1 \\
\partial p/\partial s_2
\end{bmatrix}
```

---

## `SurfaceGrad(Pu)` for vector fields

For a 3D vector field defined on a surface:

| Operator          | Size | Component ordering                     |
| ----------------- | ----:| -------------------------------------- |
| `SurfaceGrad(Pu)` | `6`  | `[ux,1, ux,2, uy,1, uy,2, uz,1, uz,2]` |

Equivalent matrix form:

```math
\nabla_\Gamma u =
\begin{bmatrix}
u_{x,1} & u_{x,2} \\
u_{y,1} & u_{y,2} \\
u_{z,1} & u_{z,2}
\end{bmatrix}
```

---

## `SurfaceSymGrad(Pu)` for vector fields

`SurfaceSymGrad(Pu)` returns the membrane strain vector
in the local tangent coordinate system.

| Operator             | Size | Component ordering |
| -------------------- | ----:| ------------------ |
| `SurfaceSymGrad(Pu)` | `3`  | `[ε11, ε22, γ12]`  |

Explicitly:

```math
\operatorname{SurfaceSymGrad}(u)
=
\begin{bmatrix}
u_{1,1} \\
u_{2,2} \\
u_{1,2} + u_{2,1}
\end{bmatrix}
```

where:

* `1` and `2` are the local tangent directions.

The operator returns membrane strains only.
No bending terms are included.

---

## `SurfaceDiv(Pu)` for vector fields

| Operator         | Size | Component ordering |
| ---------------- | ----:| ------------------ |
| `SurfaceDiv(Pu)` | `1`  | `[ux,1 + uy,2]`    |

or equivalently:

```math
\nabla_\Gamma \cdot u
```

---

# Directional / axial operators

## `AxialGrad`

Directional gradient operator along a prescribed axial direction.

This operator computes derivatives projected onto a specified axis.

Mathematically:

```math
\nabla_a u = a \cdot \nabla u
```

where:

* (a) is the prescribed axial direction.

Unlike `SurfaceGrad`, the operator does not derive its directions from the local surface geometry.

Typical applications:

* beam/spar-like formulations,
* fiber-reinforced materials,
* directional constitutive laws,
* anisotropic transport,
* projected strain operators.

---

## `AxialGrad(P)` for scalar fields

| Operator       | Size | Component ordering |
| -------------- | ----:| ------------------ |
| `AxialGrad(P)` | `1`  | `[p,a]`            |

Equivalent form:

```math
\frac{\partial p}{\partial a}
```

---

## `AxialGrad(Pu)` for vector fields

| Operator        | Size | Component ordering   |
| --------------- | ----:| -------------------- |
| `AxialGrad(Pu)` | `3`  | `[ux,a, uy,a, uz,a]` |

Equivalent form:

```math
\frac{\partial u}{\partial a}
```

## `TangentialGrad(P)` for scalar fields

Tangential derivative along a 1D embedded manifold.

| Operator            | Size | Component ordering |
| ------------------- | ----:| ------------------ |
| `TangentialGrad(P)` | `1`  | `[p,s]`            |

where:

* `s` denotes the local tangent direction.

Mathematically:

```math
\nabla_t p
=
\frac{\partial p}{\partial s}
```

---

## `TangentialGrad(Pu)` for vector fields

| Operator             | Size | Component ordering   |
| -------------------- | ----:| -------------------- |
| `TangentialGrad(Pu)` | `3`  | `[ux,s, uy,s, uz,s]` |

Equivalent form:

```math
\frac{\partial u}{\partial s}
```

## Weak-Form DSL Operators

```@docs
Grad
Div
Curl
SymGrad
ε
Id
TensorDiv
Adv
AxialGrad
TangentialGrad
SurfaceGrad
SurfaceDiv
SurfaceSymGrad
```

## Weak-Form Integration

```@docs
∫
∫Ω
∫Γ
```

## Multifield Solver

```@docs
solveField
solveEigenFields
```
