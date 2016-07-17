abstract AbstractIVP{T,Y}

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

"""

Explicit ODE representing the problem

`dy = F(t,y)` with `y(t0)=y0`

- t0, y0: initial conditions
- F!: in place version of `F` called by `F!(t,y,dy)`
- J!: (optional) computes `J=dF/dy` in place, called with `J!(t,y,J)`

"""
typealias ExplicitODE{T,Y} IVP{T,Y,Function,Void,Function}
@compat (::Type{ExplicitODE}){T,Y}(t0::T,
                                   y0::Y,
                                   F!::Function;
                                   J!::Function = forward_jacobian!(F!,similar(y0))) =
                                       ExplicitODE{T,Y}(t0,y0,similar(y0),F!,nothing,J!)


"""

Implicit ODE representing the problem

`F(t,y,dy)=0` with `y(t0)=y0` and optionally `y'(t0)=dy0`

- t0, y0: initial conditions
- F!: in place version of `F` called by `F!(t,y,dy)`
- J!: (optional) computes `J=dF/dy+a*dF/dy'` for prescribed `a`, called with `J!(t,y,dy,a)`

"""
typealias ImplicitODE{T,Y} IVP{T,Y,Void,Function,Function}
@compat (::Type{ImplicitODE}){T,Y}(t0::T,
                                   y0::Y,
                                   G!::Function;
                                   J!::Function = forward_jacobian_implicit!(F!,similar(y0)),
                                   dy0::Y = zero(y0)) =
                                       ImplicitODE{T,Y}(t0,y0,dy0,nothing,G!,J!)

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
type Step{T,S}
    t ::T
    y ::S
    dy::S
end


function show(io::IO, state::Step)
    println("t  =$(state.t)")
    println("y  =$(state.y)")
    println("dy =$(state.dy)")
end


"""

Options for ODE solvers.  This type has a key-word constructor which
will fill the structure with default values.

General:

- initstep ::T  initial step
- tstop    ::T  end integration time
- reltol   ::T  relative tolerance (m3: could this be a vector?)
- abstol   ::T  absolute tolerance (m3: could this be a vector?)
- minstep  ::T  minimal allowed step
- maxstep  ::T  maximal allowed step
- norm           function to calculate the norm in step control
- maxiters ::T  maximum number of steps
- isoutofdomain::Function checks if the solution became non-numeric (NaN or Inf)

Dense output options:

- tspan    ::Vector{T}  output times
- points   ::Symbol which points are returned: `:specified` only the
  ones in tspan or `:all` which includes also the step-points of the solver.
- stopevent  Stop integration at a zero of this function
- roottol    TODO

"""
type Options{T}
    # stepper options
    initstep ::T
    tstop    ::T
    reltol   ::T
    abstol   ::T
    minstep  ::T
    maxstep  ::T
    norm     ::Function
    maxiters ::T

    isoutofdomain::Function

    # dense output options
    tspan    ::AbstractVector{T}
    points   ::Symbol

    # m3: I think this should be an array of functions.  Depending on some
    # flag each one returns, the iteration stops or continues.  Rename it
    # to eventfns.  I like matlabs interface.
    # [value,isterminal,direction] = myEventsFcn(t,y,dy)
    # The value gets stored.
    stopevent::Function
    roottol  ::T

    function Options(;
                     tspan    = T[Inf],
                     tstop    = tspan[end],
                     reltol   = eps(T)^T(1//3)/10,
                     abstol   = eps(T)^T(1//2)/10,
                     minstep  = 10*eps(T),
                     maxstep  = 1/minstep,
                     # TODO: we need a better guess here, possibly
                     # overwrite it in the call to solve()
                     initstep = minstep,
                     norm     = Base.norm,
                     maxiters = T(Inf),
                     points   = :all,
                     stopevent = (t,y)->false,
                     roottol  = eps(T)^T(1//3),
                     isoutofdomain = isnan,
                     kargs...)
        if all(points .!= [:specified,:all])
            error("Option points = $points is not supported, use :specified or :all")
        end
        #TODO iterate over fields here?
        new(initstep,tstop,reltol,abstol,minstep,maxstep,norm,maxiters,isoutofdomain,sort(tspan),points,stopevent,roottol)
    end

end

function show{T}(io::IO, opts :: Options{T})
    for name in fieldnames(opts)
        @printf("%-20s = %s\n",name,getfield(opts,name))
    end
end


"""

This is an iterable type, each call to next(...) produces a next step
of a numerical solution to an ODE.

- ode: is the prescrived ode, along with the initial data
- stepper: the algorithm used to produce subsequent steps
- options: options passed to the stepper

"""
immutable Solver{O<:AbstractIVP,S<:AbstractStepper,T}
    ode     :: O
    stepper :: S
    options :: Options{T}
end

# filter the wrong combinations of ode and stepper
solve{T,S}(ode::T, stepper::S, options) = error("The $S doesn't support $T")


# normally we return the working array, which changes at each step and
# expect the user to copy it if necessary.  In order for collect to
# return the expected result we need to copy the output at each step.
function collect{T,S}(t::Type{Tuple{T,S}}, s::Solver)
    if maximum(s.options.tspan) == T(Inf)
        error("Attempting to collect an infinite list, use tstop or tspan with finite numbers only")
    end
    collect(t, imap(x->deepcopy(x),s))
end

function collect(s::Solver)
    if maximum(s.options.tspan) == Inf
        error("Attempting to collect an infinite list, use tstop or tspan with finite numbers only")
    end
    collect(imap(deepcopy,s))
end
