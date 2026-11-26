# Workflow

- Geometry

- FE Mesh (Gmsh)

- Fields

- Weak form

- system matrix, right hand side

- Boundary conditions

- Solver

- Resulting fields

- Post Processor

- Visualization (Gmsh)

# Setup

```julia
using LowLevelFEM
```

# Geometry

Created by Gmsh or imported STEP or IGES files.

```julia
openGeometry("part.geo")
openPreProcessor()
```

# Define the problem

```julia
mat = Material("part". E=2e5, ν=0.3)
# 3D vector field
U = Field([mat], type=:VectorField, dim=3, fieldName=:u)
# 2D scalar field
P = Field([mat], type=:ScalarField, dim=2, fieldName=:p)
```

- `type` $\to$ `:ScalarField`, `:VectorField`, `:TensorField`
- `dim` $\to$ 1,2,3
- `fieldName` $\to$ arbitrary `Symbol`

# Weak form

Solve the system matrix and right hand side

K = ∫(`<operator>` ⋅ `<matrix chain>` ⋅ `<operator>`, `Ω` or `Γ` =`<region>`, `gauss`=`<gp>`)

- `<operator>` $\to$ `U`, `Grad(U)`, `SymGrad(U)`, `Div(U)`, `Curl(U)`

- `<matrix chain>` $\to$ `Matrix` of `Number` of `ScalarField`, or a single `Number` or `ScalarField`.

- `<region>` $\to$ `String`, the name of a physical group in Gmsh. Integration domain.

- `<gp>` $\to$ `:full`, `:reduced`, `Integer`.
