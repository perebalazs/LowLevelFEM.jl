# Legacy Poisson API

```@meta
CurrentModule = LowLevelFEM
```

```@note
The Poisson operator interface is a legacy single-field API kept for compatibility.
For new models, prefer the multifield weak-form DSL (`Grad`, `Div`, `SymGrad`, `Id`, `TensorDiv`, `Adv`, `∫`) documented in [Multifield](multifield.md).
```

## Main Poisson Operators

```@docs
poissonMatrix
reactionMatrix
advectionMatrix
gradDivMatrix
symmetricGradientMatrix
curlCurlMatrix
gradMatrix
navierStokesAdvectionMatrix
tensorLaplaceMatrix
traceLaplaceMatrix
beltramiMichellMatrix
tensorDivDivMatrix
loadTensor
```

## Legacy Kernel Integrals

```@docs
∫∇oN_c_∇oN_dΩ
∫∇N_c_∇N_dΩ
∫∇xN_c_∇xN_dΩ
∫N_c_N_dΩ
∫N_c_∂N∂x_dΩ
∫N_c_∂N∂y_dΩ
∫N_c_∂N∂z_dΩ
```

`sourceVector` and `∫N_c_dΩ` are aliases of the legacy
[loadVector](linear.md) interface. Their signatures and behavior are identical.
See the [migration map](legacy.md) for scalar diffusion, reaction and advection.
