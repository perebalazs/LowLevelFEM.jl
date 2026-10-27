# Boundary Conditions and Loads

## Essential conditions

Constrain named components on physical groups:

```julia
bc = BoundaryCondition("left", field=U, ux=0.0, uy=0.0)
u = solveField(K, f; support=[bc])
```

`U` is the field descriptor used during assembly. Component names follow the
unknown name: `fieldName=:u` gives `ux,uy` for a 2D vector, while a scalar
`fieldName=:T` uses `T=...`. Values may be numbers, spatial functions or
compatible scalar fields. In mixed systems, `field=...` identifies the unknown;
a group name alone cannot disambiguate several fields on the same mesh.

## Natural loads and Robin terms

Loads in a DSL model are assembled linear forms:

```julia
f_body = ∫(U ⋅ [0.0, -1.0])
f_boundary = ∫(U ⋅ [1.0, 0.0]; Γ="right")
f = f_body + f_boundary
```

For scalar diffusion, `∫(P ⋅ q; Γ="right")` represents the boundary RHS with
`q` defined according to your weak-form sign convention. A Robin term contributes
both a boundary matrix and a load; for a 2D scalar unknown, for example:

```julia
K_robin = ∫(P ⋅ h ⋅ P; Γ="top")
f_robin = ∫(P ⋅ (h*T∞); Γ="top")
```

`Ω` and `Γ` select integration groups. Declaring an engineering `load(...)` or
`heatConvection(...)` object does not automatically add these DSL terms.

## MPCs and coordinates

Use `MPC`/`MultiPointConstraint` for master/slave relations and rigid or periodic
coupling. With local coordinate systems, constraints and ordinary RHS components
are local; use `Global(f)` for a Cartesian RHS. Returned solutions are global.

See [Boundary conditions and MPC API](../reference/boundary-conditions.md),
[periodic examples](../reference/multifield.md#Periodic-boundary-conditions), and
[Coordinate systems](coordinate-systems.md).

## Engineering declarations

`displacementConstraint`, `temperatureConstraint`, `load`, `elasticSupport`,
`heatFlux`, `heatSource`, and `heatConvection` are part of the retained
[engineering API](../reference/legacy.md). Their accepted components and effects
are documented individually.
