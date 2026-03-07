# Poisson API

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
∫N_c_dΩ
```
