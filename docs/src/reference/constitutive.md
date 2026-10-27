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

```@warning
The current Tensors extension returns tangents in `[xx,yy,zz,yz,xz,xy]` order,
while `SymGrad` and `tensorToVector` use `[xx,yy,zz,xy,yz,zx]`. For a symmetric
constitutive tensor use `p=[1,2,3,6,4,5]` and `D_dsl=D[p,p]` before combining
the tangent with those operators. Full tensor stress storage needs no such
permutation. The source implementation is unchanged.
```

## Tensor components

```@docs
tensorToMatrix
tensorToVector
```

Voigt extraction does not impose symmetry or multiply shear entries by two.
Use constitutive and kinematic conventions consistently, especially for
anisotropic materials.
