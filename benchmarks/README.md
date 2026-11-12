# Standalone verification benchmarks

Run from the root of the LowLevelFEM.jl checkout. Use the project environment so
that Julia loads LowLevelFEM from that checkout. Instantiate the project once.

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
export OPENBLAS_NUM_THREADS=1
julia --startup-file=no --threads=1 --project=. benchmarks/test-elasticity-patch.jl
julia --startup-file=no --threads=1 --project=. benchmarks/test-heat-poisson.jl
julia --startup-file=no --threads=1 --project=. benchmarks/test-nafems-le1.jl
```

Each command runs one case in a fresh Julia process. If the default Julia
environment already loads this checkout, you can also run
`julia test-nafems-le1.jl` from this directory.
All meshes are generated in memory. The scripts ignore local Gmsh configuration
and use one Gmsh thread. They require no GUI, mesh downloads or extra packages.

Each script prints the Julia, LowLevelFEM and Gmsh versions and its numerical
checks. It exits with code 0 and prints `PASS` when all checks pass. An exceeded
tolerance or non-finite error produces `FAIL` and a nonzero exit code. Runtime
errors also propagate as nonzero exit codes. The scripts run on demand,
outside `Pkg.test()` and CI.

## Analytic checks

### Plane-stress patch test

Linear triangles cover a 1 m × 1 m square of thickness 0.2 m. The analytic
displacement in metres is
`u_x = 0.001(1 + 2x + 3y)`, `u_y = 0.001(−2 + x − y)`.
With `E = 210 GPa` and `ν = 0.3`, the constant strains give all stress
components by the plane-stress Hooke law. The left edge has prescribed
displacements, and the remaining edges carry the exact tractions `σn`.

The script compares every nodal displacement and recovered stress component.
Each maximum absolute error is normalized by the largest absolute reference
component of its field. The limits are **1e-10** for displacement and
**1e-9** for stress. Linear triangles represent the affine field exactly,
so these limits allow for numerical roundoff and mesh integration.

### Steady heat conduction

On the 1 m × 1 m square, the exact temperature is
`T = sin(πx) sin(πy) K`, with `T = 0` on the boundary and thermal
conductivity `k = 3 W/(m·K)`. The volumetric source is `q = 2π²kT`,
obtained from `−kΔT = q`. The script compares temperatures at every mesh
node on 8 × 8 and 16 × 16 bilinear quadrilateral grids.

The maximum absolute nodal error on the fine grid must be **at most 0.004 K**
(0.4% of the exact peak). The fine-grid error must be **at most one third**
of the coarse-grid error, allowing for deviation from the factor of four
expected for second-order convergence when halving the element size.

## NAFEMS LE1 · Elliptic membrane

The model uses the first-quadrant section of the elliptic membrane. The inner
semi-axes are 2 m and 1 m, the outer semi-axes 3.25 m and 2.75 m, and the
thickness 0.1 m. The plane-stress material has `E = 210 GPa` and `ν = 0.3`.
Symmetry constrains `u_x` at `x = 0` and `u_y` at `y = 0`. The outer
ellipse carries a uniform outward normal traction of 10 MPa. The hole
boundary is traction-free. The published target is **σyy = 92.7 MPa**
at the inner boundary point D, `(2, 0) m`.
[NAFEMS R0006](https://www.nafems.org/publications/resource_center/r0006/) ·
[COMSOL LE1 geometry and reference value](https://doc.comsol.com/6.4/doc/com.comsol.help.models.sme.elliptic_membrane/elliptic_membrane.html)

This implementation uses its own Gmsh meshes with 12 × 8 and 24 × 16
biquadratic quadrilaterals (circumferential × radial). It compares the
recovered nodal stress at D with the published target. The relative error
`abs(σyy − 92.7 MPa) / 92.7 MPa` must be **at most 5%** on the coarse mesh
and **at most 2%** on the fine mesh. These are local acceptance limits
for this discretization.

## Recorded run

Julia 1.13.3, LowLevelFEM 1.16.2 (commit `8d732c6b`), Gmsh 4.15.2.

| Script | Observed checks | Result |
| --- | --- | --- |
| `test-elasticity-patch.jl` | Normalized max errors 5.93e-15 (displacement) and 1.84e-14 (stress) | PASS |
| `test-heat-poisson.jl` | Max nodal errors 0.0129498 K and 0.00321895 K, fine/coarse ratio 0.24857 | PASS |
| `test-nafems-le1.jl` | 92.67435 MPa (0.0277% error) and 92.82591 MPa (0.1358% error) | PASS |

The coarse LE1 mesh happens to be closer to the rounded target in this run.
These results cover the three models and the specified meshes.
