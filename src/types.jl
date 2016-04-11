abstract AbstractODE

type ExplicitODE <: AbstractODE
    t0; y0                      # initial data
    F :: Function               # solve y'=F(t,y)
    jac :: Function             # optional jacobian of F
end

type ExplicitODEInPlace <: AbstractODE
    t0; y0
    F! :: Function
    jac! :: Function
end


# plug in the numerical jacobian if none is provided
ExplicitODE(t,y,F)=ExplicitODE(t,y,F,(t,y)->fdjacobian(F,t,y))


function convert(::Type{ExplicitODEInPlace}, ode :: ExplicitODE)
    function F!(t,y,dy)
        dy[:] = ode.F(t,y)
    end
    function jac!(t,y,J)
        J[:] = ode.jac(t,y)
    end
    return ExplicitODEInPlace(ode.t0,ode.y0,F!,jac!)
end


function convert(::Type{ExplicitODE}, ode :: ExplicitODEInPlace)
    function F(t,y)
        dy = deepcopy(y)
        ode.F!(t,y,dy)
        return dy
    end
    function jac(t,y)
        n  = length(y)
        J  = Array(eltype(dy),n,n)
        ode.jac!(t,y,J)
        return J
    end
    return ExplicitODE(ode.t0,ode.y0,F,jac)
end


abstract AbstractStepper
abstract AbstractState


# this might suffice for some solvers
type Step{T,S} <: AbstractState
    t  :: T
    y  :: S
    dy :: S
end


# for debugging
function show(io::IO, state :: Step)
    println("t  =$(state.t)")
    println("y  =$(state.y)")
    println("dy =$(state.dy)")
end


# general purpose options
immutable Options{T}
    # stepper options
    initstep :: T
    tstop    :: T
    reltol   :: T
    abstol   :: T
    minstep  :: T
    maxstep  :: T
    norm     :: Function

    # dense output options
    tspan    :: Vector{T}
    points   :: Symbol
    stopevent :: Function
    roottol  :: T

    function Options(;
                     tstop    = T(Inf),
                     reltol   = eps(T)^(1/3),
                     abstol   = reltol,
                     minstep  = 10*eps(T),
                     maxstep  = 1/minstep,
                     initstep = reltol, # TODO: we need a better guess
                                        # here, possibly overwrite it
                                        # in the call to solve()
                     norm     = Base.norm,
                     tspan = [tstop],
                     points = :all,
                     stopevent = (t,y)->false,
                     roottol = eps(T)^(1/3),
                     kargs...)
        new(initstep,tstop,reltol,abstol,minstep,maxstep,norm,tspan,points,stopevent,roottol)
    end

end


function show{T}(io::IO, opts :: Options{T})
    println("")
    println("Options{$T}")
    println("tstop    = $(opts.tstop)")
    println("reltol   = $(opts.reltol)")
    println("abstol   = $(opts.abstol)")
    println("minstep  = $(opts.minstep)")
    println("maxstep  = $(opts.maxstep)")
    println("initstep = $(opts.initstep)")
    println("norm     = $(opts.norm)")
    println("tspan    = $(opts.tspan)")
    println("points   = $(opts.points)")
    println("stopevent= $(opts.stopevent)")
    println("roottol  = $(opts.roottol)")
end


# default to floating point precision
Options(args...) = Options{Float64}(args...)


# solution is a collection of an equation, an integration method
# (stepper) and its options
type Solution{T<:AbstractStepper}
    ode     :: AbstractODE
    stepper :: T
    options :: Options
end


# TODO: is this the right way to implement the mid level interface?
solve(ode, stepper; kargs...) = solve(ode, stepper, Options(kargs...))

# filter the wrong combinations of ode and stepper
solve{T,S}(ode :: T, stepper :: S, options :: Options) = error("The $S doesn't support $T")

# TODO: this might not be necessary, in the long run we should make
# ExplicitODE the in-place one and use specific constructors

# always convert ExplicitODE to ExplicitODEInPlace
solve(ode :: ExplicitODE, stepper, options :: Options) = solve(convert(ExplicitODEInPlace,ode), stepper, options)


# normally we return the working array, which changes at each step and
# expect the user to copy it if necessary.  In order for collect to
# return the expected result we need to copy the output at each step.
function collect{T}(t::Type{T}, s::Solution)
    if s.options.tstop == Inf || s.options.tspan[end] == Inf
        error("Trying to collect an infinite list")
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

function fdjacobian(F, t, x)
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
