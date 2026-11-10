# Constitutive Models and Tensor Conversion

```@meta
CurrentModule = LowLevelFEM
```

These utilities are shared by modern weak-form formulations, even though the
energy-based implementations live in `nonlinear.jl`. The dedicated nonlinear
assembly routines are [legacy](nonlinear.md).

## Isotropic linear elasticity

```@docs
ConstitutiveMatrix
D
```

The material-matrix selector is `:Axisymmetric`, while the predefined problem
type is spelled `:AxiSymmetric`. The coefficient alone does not add the
axisymmetric hoop strain or radial integration weight to a Cartesian DSL form.

## Energy-based material response

Load `Tensors` explicitly for differentiation of an energy function `ψ(C,p)`.
`C` is the right Cauchy–Green tensor. Stress and tangent are material response
fields; these functions do not run a nonlinear equilibrium iteration.

```@docs
energyMaterial
IIPiolaKirchhoff
tangentMatrix
```

Stress extraction, the energy-based tangent and `SymGrad` share the order
`[xx,yy,zz,xy,yz,zx]`. The tangent maps engineering Green-strain increments
to stress increments; use it directly without permutation or extra shear factors.

```@warning
Earlier versions of the Tensors extension used `[xx,yy,zz,yz,xz,xy]`.
Remove the old manual permutation `p=[1,2,3,6,4,5]` from existing code;
applying it again now changes the material response.
`materialTangentMatrix` also accepts the shared order on both matrix axes.
For a manually supplied matrix stored in its old `[xx,yy,zz,yz,xz,xy]` order,
convert that matrix once with `p=[1,2,3,6,4,5]; C=C[p,p]`.
```

## Tensor components

```@docs
tensorToMatrix
tensorToVector
```

Voigt extraction does not impose symmetry or multiply shear entries by two.
Use constitutive and kinematic conventions consistently, especially for
anisotropic materials.
