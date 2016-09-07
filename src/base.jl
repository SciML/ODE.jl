# The main types:
# - IVP -- holds the mathematical aspects of a IVP
# - AbstractIntegrator -- an integrator/solver  (maybe AbstractIntegrator?)
# - Problem -- holds IVP + Integrator (maybe ProblemSpec, Problem, Spec?)
# - AbstractState -- holds the iterator state
#   - Step -- holds the state at one time
# -


"""
    AbstractIVP{T,Y}

The abstract supertype of all IVPs (initial value problems).  The type
parameters `T` and `Y` correspond to the types of time and state
variable respectively.

"""
abstract AbstractIVP{T,Y}

"""
The "elements" of AbstractIVP are `t,y,dy` so the `eltype` returns `Tuple{T,Y,Y}`
"""
Base.eltype{T,Y}(::Type{AbstractIVP{T,Y}}) = Tuple{T,Y,Y}

"""

    IVP{T,Y,F,G,J} <: AbstractIVP{T,Y}

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


"""

    typealias ExplicitODE{T,Y,F,J} IVP{T,Y,F,Void,J}

Can be constructed by calling

    ODE.ExplicitODE(t0,y0,F!;J!=jacobian))

Explicit ODE representing the problem

`dy = F(t,y)` with `y(t0)=y0`

- t0, y0: initial conditions
- F!: in place version of `F` called by `F!(t,y,dy)`
- J!: (optional) computes `J=dF/dy` in place, called with `J!(t,y,J)`

"""
typealias ExplicitODE{T,Y,F,J} IVP{T,Y,F,Void,J}
@compat function (::Type{ExplicitODE}){T,Y}(t0::T,
                                            y0::Y,
                                            F!::Function;
                                            J!::Function = forward_jacobian!(F!,copy(y0)),
                                            kargs...)
    # precompute y'
    dy0 = copy(y0)
    F!(t0,y0,dy0)
    IVP(t0,y0,dy0,F!,nothing,J!)
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
typealias ImplicitODE{T,Y,G,J} IVP{T,Y,Void,G,J}
@compat function (::Type{ImplicitODE}){T,Y}(t0::T,
                                            y0::Y,
                                            G!::Function;
                                            J!::Function = forward_jacobian_implicit!(G!,similar(y0)),
                                            dy0::Y = zero(y0),
                                            kargs...)
    IVP(t0,y0,dy0,nothing,G!,J!)
end

"""

    AbstractSolver

The supertype of anything which can get you to a solution of a IVP.
Subtypes include: `AbstractIntegrator`s but also `DenseOutput`

"""
abstract AbstractSolver

Base.length(s::AbstractSolver) = error("`length` is not defined for $(typeof(s)).")

"""

The fallback generator of an abstract solver, throws an error if the
solver cannot be generated for the given initial value problem type.

"""
@compat (::Type{S}){S<:AbstractSolver}(ivp;opts...) =
    error("The solver $S doesn't support IVP of form $(typeof(ivp))")


"""

The abstract type of the actual algorithm to solve an IVP.

"""
abstract AbstractIntegrator{T} <: AbstractSolver


"""

AbstractState keeps the temporary data (state) for the iterator
Problem{::AbstractIntegrator}.

"""
abstract AbstractState{T,Y}

"""
Returns variables returned during iterations.

output(st::AbstractState) = t,y,dy
"""
output(st::AbstractState) = st.step.t, st.step.y, st.step.dy


"""

Holds a value of a function and its derivative at time t.  This is
usually used to store the solution of an IVP at particular times.

"""
type Step{T,Y}
    t ::T
    y ::Y
    dy::Y
end

output(s::Step) = s.t, s.y, s.dy

function show(io::IO, state::Step)
    println("t  =$(state.t)")
    println("y  =$(state.y)")
    println("dy =$(state.dy)")
end


"""

This is an iterable type, each call to next(...) produces a next step
of a numerical solution to an IVP.

- ivp: is the prescrived ivp, along with the initial data
- solver: the algorithm used to produce subsequent values from the ivp

"""
immutable Problem{O<:AbstractIVP,S<:AbstractSolver}
    ivp   ::O
    solver ::S
end

Base.length(prob::Problem) = length(prob.solver)

Base.eltype{O,S}(::Type{Problem{O,S}}) = eltype(O)

"""
Makes some generic operations on iterators work, like
generator comprehensions:

    tgen=(t for (t,y) in sol)
    tout=collect(tgen)

or
    errgen=(y-[exp(t)] for (t,y) in sol)
    errout=collect(errgen)
"""
Base.iteratorsize{O,S}(::Type{Problem{O,S}}) = Base.SizeUnknown()

"""
    iterate(ivp::IVP; solver=RKIntegratorAdaptive{:rk45}, opts...)

Iterate creates an iterable `Problem` instance from an `IVP` instance
(specifying the math) and from a `Type{AbstractSolver}` (the numerical
integrator).  The simplest use case is

    for (t,y,dy) in iterate(...)
        # do something with t, y an dy
    end

If the integration interval, defined by the keyword argument `tstop`,
is finite you can request all the results at once by calling

    collect(iterate(...)) # => Vector{Tuple{T,Y,Y}}

Notes:

- usually solvers require the ivp to be in a certain form, say an
 `ExplicitODE`.
- the second argument is the *Type* of the solver and not an instance.
  The instance of the solve can only be created together with the
  `ivp` as their type parameters need to match.

Input:

- `ivp::IVP`
- `S::Type{AbstractSolver}`

Output:

- `::Problem`

"""
function iterate{S<:AbstractSolver}(ivp::IVP;
                                    solver::Type{S} = RKIntegratorAdaptive{:rk45},
                                    opts...)
    Problem(ivp, solver(ivp; opts...))
end

"""
    solve(ivp::IVP; solver=RKIntegratorAdaptive{:rk45}, opts...)

Solve the initial value problem `ivp` using an algorithm `solver`
(defaults to Runge-Kutta (4,5) integrator).  One can pass additional
options to the `solver` via keyword arguments to `solve` (here denoted
as `options`).  The output is a `Solution` type (currently simply a
tuple of vectors `(Vector{T},Vector{Y})`, where `T,Y=eltype(ivp)`).

"""

function solve(ivp::IVP; opts...)

    # TODO: perhaps there is a more
    # graceful way to treat these
    # cases.  We only care about
    # infinite `tstop` but if
    # `tstop` is unspecified it
    # defaults to `tout[end]`, with
    # `tout` defaulting to
    # `[Inf]`.  Maybe we should add
    # `tout(::Problem)` to fix this?
    # Or maybe we could store
    # `tstop` in `Problem`.

    dopts = Dict(opts)
    if !in(:tstop,keys(dopts)) & !in(:tout,keys(dopts))
        error("Neither `tstop` nor `tout` was specified.")
    end
    if in(:tstop,keys(dopts)) & !isfinite(dopts[:tstop])
        error("Trying to integrate over an infinite time span, try specifying `|tstop|<Inf`.")
    end

    prob = iterate(ivp; opts...)
    T,Y = eltype(prob).parameters
    tout = Array(T,0)
    yout = Array(Y,0)
    dyout = Array(Y,0)
    for (t,y,dy) in prob
        push!(tout,t)
        push!(yout,copy(y))
        push!(dyout,copy(dy))
    end
    return Solution{T,Y}(tout,yout,dyout)
end


"""

Stores a solution to the `ivp`.

"""
immutable Solution{T,Y}
    t::Vector{T}
    y::Vector{Y}
    dy::Vector{Y}
end

"""
Experimental support for the interpolation
"""
function interpolate(sol::Solution, t)
    if (t>sol.t[end]) | (t<sol.t[1])
        error("Out of range")
    else
        i = findfirst(ti -> t>ti, sol.t)
        theta = (t-sol.t[i])/(sol.t[i+1]-sol.t[i])
        return sol.y[i]+theta*(sol.y[i+1]-sol.y[i])
    end
end


# In Julia 0.5 the collect needs length to be defined, we cannot do
# that for a Problem but we can implement our own collect
function collect(prob::Problem)
    pairs = Array(eltype(prob),0)
    for (t,y,dy) in prob
        push!(pairs,(t,copy(y),copy(dy)))
    end
    return pairs
end


# Iteration: take one step on a IVP `Problem`
#
# Defines:
# start(iter) -> state
# next(iter, state) -> output(state), state
# done(iter, state) -> bool
#
# Perhaps unintuitively, the next step is computed in `done`.  Such
# implementation allows to decide if the iterator is exhausted in case
# when the next step was computed but it was deemed incorrect.  In
# such situation `done` returns `false` after computing the step and
# the failed step never sees the light of the day (by not being
# returned by `next`).
#
# TODO: this implementation fails to return the zeroth step (t0,y0)
#
# TODO: store the current Step outside of the actual state

Base.start(prob::Problem) = init(prob.ivp, prob.solver)

function Base.done(prob::Problem, st)
    # Determine whether the next step can be made by calling the
    # stepping routine.  onestep! will take the step in-place.
    status = onestep!(prob.ivp, prob.solver, st)
    if status==cont
        return false
    elseif status==finish
        return true
    else #if status==abort
        warn("aborting")
        return true
    # else
    #     error("unsported Status: $status")
    end
end

function Base.next(prob::Problem, st)
    # Output the step (we know that `done` allowed it, so we are safe
    # to do it)
    return output(st), st
end

"""

Holds the solver status, used inside of `onestep!`.

Values:

- cont -- continue integration
- abort -- abort integration
- finish -- integration reached the end

Statuses can be combined with &:
- cont&cont == cont
- finish&cont == finish
- abort&cont == abort
- abort&finish = abort

"""
@enum Status cont=1 abort=0 finish=-1
# The values of Statuses are chose to turn & into a *:
@compat Base.:&(s1::Status, s2::Status) = Status(Int(s1)*Int(s2))

#####
# Interface to implement by solvers to hook into iteration
#####
#
# See runge_kutta.jl and rosenbrock.jl for example implementations.

# A stepper has to implement
# - init
# - output
# and either
# - onestep!
# - trialstep!, errorcontrol! and accept!

"""

Take a step, modifies `state` in-place.  This is the core function to
be implemented by a solver.  However, if possible solvers should opt
to implement the sub-step functions `trialstep!`, `errorcontrol!` and
`accept!`, instead of directly `onestep!`.

Input:

- prob::Problem, state::AbstractState

Output:

- Status

substeps.

"""
function onestep!(ivp::IVP, integ::AbstractIntegrator, state::AbstractState)
    opt = integ.opts
    while true
        status = trialstep!(ivp, integ, state)
        err, status_err = errorcontrol!(ivp, integ, state)
        status &= status_err
        if err<=1
            # a successful step
            status &= accept!(ivp, integ, state)
            return status
        elseif status==abort || status==finish
            return status
        end
        # if we get here: try step again with updated state (step
        # size, order).
    end
end


"""

Advances the solution by trying to compute a single step.  The new
step is kept in the `state` in work arrays so that `errorcontrol!` can
compute the magnitude of its error.  If the error is small enough
`accept!` updates `state` to reflect the state at the new time.

Returns `Status`.

"""
trialstep!{I<:AbstractIntegrator}(::IVP, ::I, ::AbstractState) =
    error("Function `trialstep!` and companions (or alternatively `onestep!`) need to be implemented for adaptive integrator $I")

"""

Estimates the error (such that a step is accepted if err<=1).
Depending on the stepper it may update the state, e.g. by computing a
new dt or a new order (but not by computing a new solution!).

Returns `(err,Status)`.

If the `status==abort` then the integration is aborted, status values
of `cont` and `finish` are ignored.

"""
errorcontrol!{I<:AbstractIntegrator}(::IVP, ::I, ::AbstractState) =
    error("Function `errorcontrol!` and companions (or alternatively `onestep!`) need to be implemented for adaptive integrator $I")

"""

Accepts (in-place) the computed step.  Called if `errorcontrol!` gave
a small enough error.

Returns `Status`.

"""
accept!{I<:AbstractIntegrator}(::IVP, ::I, ::AbstractState) =
    error("Function `accept!` and companions (or alternatively `onestep!`) need to be implemented for adaptive integrator $I")
