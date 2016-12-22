abstract ODEjlAlgorithm <: AbstractODEAlgorithm
# Making the ODE-solver functions into types lets us dispatch on them.
#  Used in the DiffEqBase interface.
immutable ode23 <: ODEjlAlgorithm end
immutable ode45 <: ODEjlAlgorithm end
immutable ode23s <: ODEjlAlgorithm end
immutable ode78 <: ODEjlAlgorithm end
immutable ode4 <: ODEjlAlgorithm end
immutable ode4ms <: ODEjlAlgorithm end
immutable ode4s <: ODEjlAlgorithm end
