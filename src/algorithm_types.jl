abstract type ODEjlAlgorithm <: AbstractODEAlgorithm end
# Making the ODE-solver functions into types lets us dispatch on them.
#  Used in the DiffEqBase interface.
struct ode23 <: ODEjlAlgorithm end
struct ode45 <: ODEjlAlgorithm end
struct ode23s <: ODEjlAlgorithm end
struct ode78 <: ODEjlAlgorithm end
struct ode4 <: ODEjlAlgorithm end
struct ode4ms <: ODEjlAlgorithm end
struct ode4s <: ODEjlAlgorithm end
