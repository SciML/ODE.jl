```@meta
CurrentModule = ODE
```

# Base

The file `base.jl` implements the most basic iterator infrastructure
for solvers and the definitions of the types representing general IVP
(initial value problem) and solvers.

## Predefined types of initial value problems

```@docs
AbstractIVP
IVP
ExplicitODE
ImplicitODE
```
## Solver architecture

```@docs
AbstractSolver
AbstractIntegrator
AbstractState
```

The fallback constructor for `AbstractSolver(ivp::IVP;opts...)` ensures
that an error is thrown if a solver is constructed for an unsupported
type of the given IVP.

## Fallback functions for solvers

```@docs
Base.length(::AbstractSolver)
output(::AbstractState)
Base.eltype{T,Y}(::Type{AbstractIVP{T,Y}})
```
