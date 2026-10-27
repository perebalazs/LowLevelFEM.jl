# Troubleshooting

## Missing regions or wrong dimensions

Check that every material, `Ω`, `Γ` and support group exists in the current
Gmsh model. Structured helpers document their generated names. Do not reuse a
field descriptor after replacing its model/mesh; reconstruct descriptors and
matrix patterns for the new mesh.

## Singular or indefinite systems

Check essential constraints, rigid modes, disconnected regions and a pressure
reference in mixed incompressible problems. `Symmetric(K)` declares symmetry;
it does not make an indefinite matrix positive definite. Use a suitable direct
solver for saddle-point systems before considering specialized iteration.

## Field shape or coefficient errors

Inspect `P.dim`, `P.pdim`, `isNodal(f)`, `isElementwise(f)` and column counts.
A 2D tensor field still stores nine components. Match element supports and
component ordering before combining fields. Some assembly paths require nodal
scalar coefficients; convert explicitly where required by the method.

## Local loads and constraints

With `coordSys`, ordinary RHS components and prescribed values are local.
Wrap Cartesian loads in `Global(f)`. In mixed systems bind each essential
condition with `field=P`. See [coordinate systems](coordinate-systems.md).

## Missing optional solver/material packages

Load `IterativeSolvers` for CG/GMRES; automatic incomplete Cholesky/ILU requires
`Preconditioners`/`IncompleteLU`. ARPACK eigen solves require `Arpack` and
energy differentiation requires `Tensors`. No optional package is loaded
implicitly. Existing installed extensions become active after `using` them.

## Eigenvalue selection

`shift` is a one-sided real lower bound; ARPACK `sigma` is a shift-invert target
and may select modes from both sides. They cannot be combined. Use `output=:raw`
for complete complex eigenpairs, or `part=:both` for real/imaginary fields.
Do not interpret every eigenvalue as a vibration frequency.

## Known interface inconsistencies

See the [public API audit](../reference/api-audit.md) for exported undefined
names, `Adv(P)`, and tangent component-order compatibility. These are documented
limitations; this documentation update does not change executable source.
