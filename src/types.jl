abstract AbstractODE


type ExplicitODE <: AbstractODE
    t0; y0
    F! :: Function
    jac! :: Function
end


function explicit_ineff(t0,y0,F;jac = (t,y)->fdjacobian(F,t,y))
    function F!(t,y,dy)
        # this is why we can't handle a scalar type any more
        dy[:] = F(t,y)
    end
    function jac!(t,y,J)
        J[:] = jac(t,y)
    end
    return ExplicitODE(t0,y0,F!,jac!)
end


abstract AbstractStepper
abstract AbstractState


# this might suffice for some solvers
type Step{T,S}
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
    maxiters :: T

    # dense output options
    tspan    :: Vector{T}
    points   :: Symbol
    stopevent :: Function
    roottol  :: T

    function Options(;
                     tstop    = T(Inf),
                     tspan = [tstop],
                     reltol   = eps(T)^(1/3)/10,
                     abstol   = eps(T)^(1/2)/10,
                     minstep  = 10*eps(T),
                     maxstep  = 1/minstep,
                     # TODO: we need a better guess here, possibly
                     # overwrite it in the call to solve()
                     initstep = max(min(reltol,abstol,maxstep),minstep),
                     norm     = Base.norm,
                     maxiters = T(Inf),
                     points = :all,
                     stopevent = (t,y)->false,
                     roottol = eps(T)^(1/3),
                     kargs...)
        if all(points .!= [:specified,:all])
            error("Option points = $points is not supported, use :specified or :all")
        end
        new(initstep,tstop,reltol,abstol,minstep,maxstep,norm,maxiters,sort(tspan),points,stopevent,roottol)
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


# normally we return the working array, which changes at each step and
# expect the user to copy it if necessary.  In order for collect to
# return the expected result we need to copy the output at each step.
function collect{T}(t::Type{T}, s::Solution)
    if any(s.options.tspan .== Inf)
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
