abstract AbstractODE

"""
Explicitly defined ODE of form dy = F(t,y).

Fields:

- t0, y0: initial conditions
- F!: ODE function `F!(t,y,dy)` which modifies `dy` in-place
- jac!: TODO
"""
immutable ExplicitODE{T,S} <: AbstractODE
    t0  ::T
    y0  ::S
    F!  ::Function
    jac!::Function
end

# TODO: change (t,y,J)->fdjacobian(F!,t,y,J) to fdjacobian!(F!)
ExplicitODE{T,S}(t0::T, y0::S, F!::Function;
                 jac!::Function = (t,y,J)->fdjacobian!(F!,t,y,J)) =
    ExplicitODE{T,S}(t0,y0,F!,jac!)


"""
This type is not yet implemented, but will serve as an implicitly
defined ODE (i.e. ODE of the form F(t,y,y')=0.
"""
immutable ImplicitODE{T,S} <: AbstractODE
end


"""
Convert a out-of-place explicitly defined ODE function to an in-place function.

Note, this does not help with memory allocations.
"""
function explicit_ineff(t0, y0, F::Function;
                        jac = (t,y)->fdjacobian(F,t,y))
    function F!(t,y,dy)
        # this is why we can't handle a scalar type any more
        copy!(dy,F(t,y))
    end
    function jac!(t,y,J)
        copy!(J,jac(t,y))
    end
    return ExplicitODE(t0,y0,F!,jac!)
end


"""
The abstract type of the actual algorithm to solve an ODE.
"""
abstract AbstractStepper


"""
AbstractState keeps the temporary data (state) for the iterator
Solver{::AbstractStepper}.
"""
abstract AbstractState

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
immutable Options{T} # m3:  T->Et
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
    tspan    ::Vector{T}
    points   ::Symbol

    # m3: I think this should be an array of functions.  Depending on some
    # flag each one returns, the iteration stops or continues.  Rename it
    # to eventfns.  I like matlabs interface.
    # [value,isterminal,direction] = myEventsFcn(t,y,dy)
    # The value gets stored.
    stopevent::Function
    roottol  ::T

    function Options(;
                     tstop    = T(Inf),
                     tspan    = T[tstop],
                     reltol   = eps(T)^T(1//3)/10,
                     abstol   = eps(T)^T(1//2)/10,
                     minstep  = 10*eps(T),
                     maxstep  = 1/minstep,
                     # TODO: we need a better guess here, possibly
                     # overwrite it in the call to solve()
                     initstep = max(min(reltol,abstol,maxstep),minstep),
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
type Solver{T<:AbstractStepper}
    ode     :: AbstractODE
    stepper :: T
    options :: Options
end


# TODO: is this the right way to implement the mid level interface?
solve(ode, stepper; kargs...) = solve(ode, stepper, Options(kargs...))

# filter the wrong combinations of ode and stepper
solve{T,S}(ode :: T, stepper :: S, options :: Options) = error("The $S doesn't support $T")


# normally we return the working array, which changes at each step and
# expect the user to copy it if necessary.  In order for collect to
# return the expected result we need to copy the output at each step.
function collect{T,S}(t::Type{Tuple{T,S}}, s::Solver)
    if maximum(s.options.tspan) == T(Inf)
        error("Attempting to collect an infinite list, use tstop or tspan with finite numbers only")
    end
    collect(t, imap(x->deepcopy(x),s))
end


# some leftovers from the previous implementation

# FIXME: This doesn't really work if x is anything but a Vector or a scalar
function fdjacobian(F, t, x::Number)
    ftx = F(t, x)

    # The 100 below is heuristic
    dx = (x .+ (x==0))./100
    dFdx = (F(t,x+dx)-ftx)./dx

    return dFdx
end

function fdjacobian(F, t, x::Vector)
    ftx = F(t, x)
    lx = max(length(x),1)
    dFdx = zeros(eltype(x), lx, lx)
    for j = 1:lx
        # The 100 below is heuristic
        dx = zeros(eltype(x), lx)
        dx[j] = (x[j] .+ (x[j]==0))./100
        dFdx[:,j] = (F(t,x+dx)-ftx)./dx[j]
    end
    return dFdx
end

function fdjacobian!{T}(F!, t, x::Vector{T}, J::Array{T,2})
    ftx = similar(x)
    ftx2= similar(x)
    dx  = similar(x)
    F!(t,x,ftx)
    lx = max(length(x),1)
    dFdx = zeros(eltype(x), lx, lx)
    for j = 1:lx
        # The 100 below is heuristic
        dx[j] = (x[j] .+ (x[j]==0))./100
        F!(t,x+dx,ftx2)
        J[:,j] = (ftx2-ftx)./dx[j]
    end
end
