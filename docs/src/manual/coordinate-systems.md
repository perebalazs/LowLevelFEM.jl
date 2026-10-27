# Coordinate Systems

Local coordinates define the component basis of the constrained solver system;
they do not move mesh nodes or change the geometry of the weak form.

In 2D, give one direction field. For an existing vector descriptor `U`:

```julia
e1 = VectorField(U, "body", [0.0, 1.0, 0.0])
cs = CoordinateSystem(e1)
f_cartesian = ∫(U ⋅ [1.0, 0.0])
bc = BoundaryCondition("left", field=U, ux=0.0, uy=0.0)
u = solveField(K, Global(f_cartesian); support=[bc], coordSys=[cs])
u_local = rotateResult(u; toCoordSys=[cs])
```

`K` is assembled in Cartesian coordinates. With `coordSys`, an unwrapped RHS
would be interpreted as already local; `Global(...)` requests the Cartesian-to-local
transformation. Prescribed component values are local. The returned `u` is in
Cartesian coordinates; `rotateResult` changes the representation for inspection.

In 3D, supply two direction fields. They must belong to the same descriptor;
the second sets the local 1–2 plane and is orthogonalized against the first.
`handedness=:right` or `:left` controls the generated remaining direction.
Direction fields may vary in space, but must be nondegenerate and have one column.
For example a radial basis is undefined at the cylinder axis.

MPCs relate local solver components; different master/slave bases can represent
rotated periodic relations. Local coordinates combined with `reducedOrder=true`
on the same field are currently unsupported; using them on different fields of
a mixed system is permitted. Check additional spectral/MPC restrictions in
[Solution API](../reference/solution.md).

See [CoordinateSystem, rotateNodes and rotateResult](../reference/boundary-conditions.md)
for exact signatures and supported field representations.
