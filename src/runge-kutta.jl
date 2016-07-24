# This file contains the implementation of explicit Runkge-Kutta
# solver from (Hairer & Wanner 1992 p.134, p.165-169).

include("tableaus.jl")

# intermediate level interface


"""

A general Runge-Kutta stepper (it can represent either, a fixed step
or an adaptive step algorithm).

"""
immutable RKStepper{Kind,Name,T,O<:Options} <: AbstractStepper{T}
    tableau::TableauRKExplicit{T}
    options::O
end


typealias RKStepperFixed    RKStepper{:fixed}
typealias RKStepperAdaptive RKStepper{:adaptive}


@compat function (::Type{RKStepper{Kind,Name,T}}){Kind,Name,T}(;options...)
    tab = convert(TableauRKExplicit{T},tableaus_rk_explicit[Name])
    if Kind == :fixed
        opts = FixedOptions{T}(;options...)
        if isadaptive(tab)
            error("Cannot construct a fixed step method from an adaptive step tableau")
        end
    elseif Kind == :adaptive
        opts = AdaptiveOptions{T}(;options...)
        if !isadaptive(tab)
            error("Cannot construct an adaptive step method from an fixed step tableau")
        end
    end
    RKStepper{Kind,Name,T,typeof(opts)}(tab,opts)
end


order(stepper::RKStepper) = minimum(order(stepper.tableau))

name(stepper::RKStepper) = stepper.tableau.name

solve{T,S<:RKStepper}(ode::ExplicitODE{T}, stepper::Type{S}; options...) =
    Solver(ode,stepper{T}(;options...))

# lower level interface

# explicit RK stepper

"""

Pre allocated arrays to store temporary data.  Used only by
Runge-Kutta stepper.

"""
type RKWorkArrays{Y}
    y   ::Y
    ynew::Y
    yerr::Y
    ks  ::Vector{Y}
end


"""
State for the Runge-Kutta stepper.
"""
type RKState{T,Y} <: AbstractState{T,Y}
    step    ::Step{T,Y}
    dt      ::T
    newdt   ::T
    work    ::RKWorkArrays{Y}
    timeout ::Int
    # This is not currently incremented with each step
    iters   ::Int
end


output(st::RKState) = st.step.t, st.step.y


function show(io::IO, state::RKState)
    show(state.step)
    println("dt      = $(state.dt)")
    println("timeout = $(state.timeout)")
    println("work    = $(state.work)")
end


function init{O<:ExplicitODE,S<:RKStepper}(s::Solver{O,S})
    stepper = s.stepper
    t0, dt0, y0 = s.ode.t0, stepper.options.initstep, s.ode.y0

    # clip the dt0 if t0+dt0 exceeds tstop
    dt0 = min(dt0,stepper.options.tstop-t0)

    lk = lengthks(s.stepper.tableau)
    work = RKWorkArrays(zero(y0), # y
                        zero(y0), # ynew
                        zero(y0), # yerr
                        Array(typeof(y0), lk)) # ks

    # we have to allocate each component separately
    for i = 1:lk
        work.ks[i]=zero(y0)
    end

    # pre-initialize work.ks[1]
    s.ode.F!(t0,y0,work.ks[1])

    step = Step(t0,copy(y0),copy(work.ks[1]))

    timeout = 0 # for step control
    return RKState(step,dt0,dt0,work,timeout,0)
end


#####################
# Fixed step method #
#####################


function onestep!{O<:ExplicitODE,S<:RKStepperFixed}(s::Solver{O,S}, state::RKState)
    step = state.step
    work = state.work

    if step.t >= s.stepper.options.tstop
        # nothing left to integrate
        return StatusFinished
    end

    dof  = length(step.y)
    b    = s.stepper.tableau.b
    dt   = min(state.dt,s.stepper.options.tstop-step.t)

    copy!(work.ynew,step.y)

    for k=1:length(b)
        calc_next_k!(work, k, s.ode, s.stepper.tableau, step, dt)
        for d=1:dof
            work.ynew[d] += dt * b[k]*work.ks[k][d]
        end
    end
    step.t += dt
    copy!(step.y,work.ynew)
    return StatusContinue
end


########################
# Adaptive step method #
########################


const timeout_const = 5

# `trialstep!` ends with a step computed for the stepsize `state.dt`
# and stores it in `work.y`, so `work.y` contains a candidate for
# `y(t+dt)` with `dt=state.dt`.
function trialstep!{O<:ExplicitODE,S<:RKStepperAdaptive}(sol::Solver{O,S}, state::RKState)
    work    = state.work
    step    = state.step
    stepper = sol.stepper
    tableau = stepper.tableau
    options = stepper.options

    # use the proposed step size to perform the computations
    state.dt = state.newdt
    dt = state.dt

    if step.t >= options.tstop
        # nothing left to integrate
        return StatusFinished
    end

    if dt < options.minstep
        # minimum step size reached
        return StatusFailed
    end

    # work.y and work.yerr and work.ks are updated after this step
    rk_embedded_step!(work, sol.ode, tableau, step, dt)

    return StatusContinue
end

# computes the error for the candidate solution `y(t+dt)` with
# `dt=state.dt` and proposes a new time step
function errorcontrol!{O<:ExplicitODE,S<:RKStepperAdaptive}(sol::Solver{O,S}, state::RKState)
    work = state.work
    step = state.step
    stepper = sol.stepper
    tableau = stepper.tableau
    timeout = state.timeout
    options = stepper.options
    err, state.newdt, state.timeout =
        stepsize_hw92!(work, step, tableau, state.dt, state.timeout, options)

    # trim in case newdt > dt
    state.newdt = min(state.newdt, options.tstop-state.step.t)

    if err > 1
        # The error is too large, the step will be rejected.  We reset
        # the timeout and set the new stepsize
        state.timeout = timeout_const
    end

    return err, StatusContinue
end

# Here we assume that trialstep! and errorcontrol! have already been
# called, that is `work.y` holds `y(t+dt)` with `dt=state.dt`, and
# error was small enough for us to keep `y(t+dt)` as the next step.
function accept!{O<:ExplicitODE,S<:RKStepperAdaptive}(sol::Solver{O,S}, state::RKState)
    work    = state.work
    step    = state.step
    tableau = sol.stepper.tableau

    # preload ks[1] for the next step
    if tableau.isFSAL
        copy!(work.ks[1],work.ks[end])
    else
        sol.ode.F!(step.t+state.dt, work.ynew, work.ks[1])
    end

    # Swap bindings of y and ytrial, avoids one copy
    step.y, work.ynew = work.ynew, step.y
    # state.dt holds the size of the last successful step
    step.t += state.dt

    return StatusContinue
end


##########################
# Lower level algorithms #
##########################

function rk_embedded_step!(work      ::RKWorkArrays,
                           ode       ::ExplicitODE,
                           tableau   ::Tableau,
                           last_step ::Step,
                           dt)
    # Does one embedded R-K step updating work.ynew, work.yerr and work.ks.
    # Assumes that work.ks[:,1] is already calculated!
    # Modifies work.y, work.ynew and work.yerr only

    y      = last_step.y
    dof    = length(y)
    b      = tableau.b

    fill!(work.ynew, zero(eltype(y)))
    fill!(work.yerr, zero(eltype(y)))

    for s=1:lengthks(tableau)
        # we skip the first step beacause we assume that work.ks[1] is
        # already computed
        if s > 1
            calc_next_k!(work, s, ode, tableau, last_step, dt)
        end
        for d=1:dof
            work.ynew[d] += b[1,s]*work.ks[s][d]
            work.yerr[d] += b[2,s]*work.ks[s][d]
        end
    end

    for d=1:dof
        work.yerr[d] = dt*(work.ynew[d]-work.yerr[d])
        work.ynew[d] = y[d] + dt*work.ynew[d]
    end

    return nothing
end


function stepsize_hw92!{T}(work,
                           last_step ::Step,
                           tableau   ::TableauRKExplicit,
                           dt        ::T,
                           timeout   ::Int,
                           options   ::Options{T})
    # Estimates the error and a new step size following Hairer &
    # Wanner 1992, p167 (with some modifications)
    #
    # If timeout>0 no step size increase is allowed, timeout is
    # decremented in here.
    #
    # Returns the error, newdt and the number of timeout-steps
    #
    # TODO:
    # - allow component-wise reltol and abstol?
    # - allow other norms

    ord = minimum(order(tableau))
    timout_after_nan = 5
    # fac = T[0.8, 0.9, (0.25)^(1/(ord+1)), (0.38)^(1/(ord+1))][1]
    fac = T(8//10)
    facmax = T(5) # maximal step size increase. 1.5-5
    facmin = 1./facmax  # maximal step size decrease. ?
    dof = length(last_step.y)

    # in-place calculate yerr./tol
    for d=1:dof

        # TODO: for some reason calling options.isoutofdomain
        # generates a lot of allocations

        # if options.isoutofdomain(work.y[d])::Bool
        if isnan(work.y[d])
            return T(10), dt*facmin, timout_after_nan
        end

        y0 = last_step.y[d] # TODO: is this supposed to be the last successful step?
        y1 = work.ynew[d]    # the approximation to the next step
        sci = (options.abstol + options.reltol*max(norm(y0),norm(y1)))
        work.yerr[d] ./= sci # Eq 4.10
    end

    # TOOD: should we use options.norm here as well?
    err   = options.norm(work.yerr) # Eq. 4.11
    newdt = min(options.maxstep, dt*clamp(fac*(1/err)^(1/(ord+1)),facmin,facmax)) # Eq 4.13 modified

    if timeout > 0
        newdt = min(newdt, dt)
        timeout -= 1
    end

    return err, newdt, timeout
end


# For clarity we pass the RKWorkArrays part of the state separately,
# this is the only part of state that can be changed here
function calc_next_k!(work      ::RKWorkArrays,
                      i         ::Int,
                      ode       ::ExplicitODE,
                      tableau   ::Tableau,
                      last_step ::Step,
                      dt)
    dof = length(last_step.y)
    t, a, c = last_step.t, tableau.a, tableau.c

    copy!(work.y,last_step.y)
    for j=1:i-1
        for d=1:dof
            work.y[d] += dt * work.ks[j][d] * a[i,j]
        end
    end
    ode.F!(t + c[i]*dt, work.y, work.ks[i])
    return nothing
end
