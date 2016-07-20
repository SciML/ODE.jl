# The main types:
# - IVP -- holds the mathematical aspects of a IVP
# - AbstractStepper -- an integrator/solver  (maybe AbstractIntegrator?)
# - Solver -- holds IVP + Stepper (maybe ProblemSpec, Problem, Spec?)
# - AbstractState -- holds the iterator state
#   - Step -- holds the state at one time
# -


abstract AbstractIVP{T,Y}
Base.eltype{T,Y}(::Type{AbstractIVP{T,Y}}) = T,Y

"""

Defines the mathematical part of an IVP (initial value problem)
specified in the general form:

`F(t, y) =  G(t, y, dy)` with `y(t0)= y0`

Depending on the combination of the parameters this type can represent
a wide range of problems, including ODE, DAE and IMEX.  Nevertheless
not all solvers will support any combinations of `F` and `G`.  Note
that not specifying `G` amounts to `G=dy/dt`.


- `tspan` -- tuple `(start_t,end_t)`
- `y0` -- initial condition
- `F!` -- in-place `F` function `F!(t,y,res)`.  If `F=0` set to `nothing`.
- `G!` -- in-place `G` function `G!(t,y,dy,res)`.  If `G=dy/dt` then
          set to `nothing` (or `dy` if the solver supports this).  Can
          also be a mass matrix for a RHS `M dy/dt`
- `J!` -- in-place Jacobian function `J!(t,y,dy,res)`.

TODO: how to fit the sparsity pattern in J?

"""
type IVP{T,Y,F,G,J} <: AbstractIVP{T,Y}
    t0  ::T
    y0  ::Y
    dy0 ::Y
    F!  ::F
    G!  ::G
    J!  ::J
end
@compat Base.eltype(t::Type{IVP}) = eltype(supertype(t))
Base.eltype(t::IVP) = eltype(typeof(t))


"""

Explicit ODE representing the problem

`dy = F(t,y)` with `y(t0)=y0`

- t0, y0: initial conditions
- F!: in place version of `F` called by `F!(t,y,dy)`
- J!: (optional) computes `J=dF/dy` in place, called with `J!(t,y,J)`

"""
typealias ExplicitODE{T,Y} IVP{T,Y,Function,Void,Function}
@compat function (::Type{ExplicitODE}){T,Y}(t0::T,
                                            y0::Y,
                                            F!::Function;
                                            J!::Function = forward_jacobian!(F!,similar(y0)))
    ExplicitODE{T,Y}(t0,y0,similar(y0),F!,nothing,J!)
end

"""

Implicit ODE representing the problem

`G(t,y,dy)=0` with `y(t0)=y0` and optionally `y'(t0)=dy0`

- t0, y0: initial conditions
- G!: in place version of `G` called by `G!(res,t,y,dy)`,
      returns residual in-place in `res`.
- J!: (optional) computes `J=dF/dy+a*dF/dy'` for prescribed `a`, called with `J!(out,t,y,dy,a)`.
      Returns Jacobian in-place in `out`.

"""
typealias ImplicitODE{T,Y} IVP{T,Y,Void,Function,Function}
@compat function (::Type{ImplicitODE}){T,Y}(t0::T,
                                            y0::Y,
                                            G!::Function;
                                            J!::Function = forward_jacobian_implicit!(G!,similar(y0)),
                                            dy0::Y = zero(y0))
    ImplicitODE{T,Y}(t0,y0,dy0,nothing,G!,J!)
end

"""

The abstract type of the actual algorithm to solve an ODE.

"""
abstract AbstractStepper{T}


"""

AbstractState keeps the temporary data (state) for the iterator
Solver{::AbstractStepper}.

"""
abstract AbstractState{T,Y}

# m3:
# - docs
# - maybe use the typevars as defined in make_consistent_types for t,
#   y, dy?  T->Et, S->Ty
#   (or something else consistent throughout, maybe nicer would be all
#   uppercase: ET, EFY, TT, TY).
# - if find `Step` a bit confusing name, in particular combined with
#   AbstractStepper, but not sure what's better.

"""

Holds a value of a function and its derivative at time t.  This is
usually used to store the solution of an ODE at particular times.

"""
type Step{T,Y}
    t ::T
    y ::Y
    dy::Y
end


function show(io::IO, state::Step)
    println("t  =$(state.t)")
    println("y  =$(state.y)")
    println("dy =$(state.dy)")
end


"""

This is an iterable type, each call to next(...) produces a next step
of a numerical solution to an ODE.

- ode: is the prescrived ode, along with the initial data
- stepper: the algorithm used to produce subsequent steps


"""
immutable Solver{O<:AbstractIVP,S<:AbstractStepper}
    ode     ::O
    stepper ::S
end
#m3:
# - calling this `Solver` still trips me up

Base.eltype{O}(::Type{Solver{O}}) = eltype(O)
Base.eltype{O}(::Solver{O}) = eltype(O)

# filter the wrong combinations of ode and stepper
solve{O,S}(ode::O, stepper::Type{S}, options...) =
    error("The $S doesn't support $O")

# In Julia 0.5 the collect needs length to be defined, we cannot do
# that for a solver but we can implement our own collect
function collect(s::Solver)
    T,Y = eltype(s)
    pairs = Array(Tuple{T,Y},0)
    for (t,y) in s
        push!(pairs,(t,copy(y)))
    end
    return pairs
end
