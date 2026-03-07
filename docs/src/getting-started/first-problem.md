# First Problem

This is the shortest end-to-end linear elasticity workflow.

```julia
using LowLevelFEM

gmsh.initialize()
gmsh.open("your_model.geo")

mat = Material("body", E=2e5, ν=0.3)
prob = Problem([mat], type=:PlaneStress)

bc = displacementConstraint("supp", ux=0, uy=0)
ld = load("load", fy=-1)

u = solveDisplacement(prob, load=[ld], support=[bc])
S = solveStress(u)

showDoFResults(u, :ux)
showStressResults(S, :s)

gmsh.finalize()
```

## Next Steps

- Continue with [Mesh and Physical Groups](mesh-and-physical-groups.md)
- See full worked models in [Tutorials](../tutorials/index.md)
