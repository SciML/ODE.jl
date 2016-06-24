# This file contains the implementation of explicit Runkge-Kutta
# solver from (Hairer & Wanner 1992 p.134, p.165-169).

include("tableaus.jl")

# intermediate level interface

"""

A general Runge-Kutta stepper (it cen represent either, a fixed step
or an adaptive step algorithm).

"""
immutable RKStepper{Kind,tab_name,T<:Number} <: AbstractStepper
    tableau::Tableau
    function RKStepper()
        tab = eval(tab_name)
        if Kind == :fixed && isadaptive(tab)
            error("Cannot construct a fixed step method from an adaptive step tableau")
        elseif Kind == :adaptive && !isadaptive(tab)
            error("Cannot construct an adaptive step method from an fixed step tableau")
        end
        new(convert(T,tab))
    end
end


typealias RKStepperFixed    RKStepper{:fixed}
typealias RKStepperAdaptive RKStepper{:adaptive}


order(stepper::RKStepper) = minimum(order(stepper.tableau))

name(stepper::RKStepper) = typeof(stepper.tableau)

# TODO: possibly handle the initial stepsize and the tableau conversion here?
solve{K,S,T}(ode::ExplicitODE, stepper::RKStepper{K,S,T}, options::Options{T}) =
    Solver{RKStepper{K,S,T}}(ode,stepper,options)


# lower level interface

# explicit RK stepper

"""

Pre allocated arrays to store temporary data.  Used only by
Runge-Kutta stepper.

"""
type RKWorkArrays{T}
    y   ::T
    ynew::T
    yerr::T
    ks  ::Vector{T}
end


"""
State for the Runge-Kutta stepper.
"""
type RKState{T,S} <: AbstractState
    step    ::Step{T,S}
    dt      ::T
    work    ::RKWorkArrays{S}
    timeout ::Int
    # This is not currently incremented with each step
    iters   ::Int
end

function show(io::IO, state::RKState)
    show(state.step)
    println("dt      = $(state.dt)")
    println("timeout = $(state.timeout)")
    println("work    = $(state.work)")
end


function start{T<:RKStepper}(s::Solver{T})
    t0, dt0, y0 = s.ode.t0, s.options.initstep, s.ode.y0

    # TODO: we should do the Butcher table conversion somewhere
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

    step = Step(t0,deepcopy(y0),deepcopy(work.ks[1]))

    timeout = 0 # for step control
    return RKState(step,dt0,work,timeout,0)
end


#####################
# Fixed step method #
#####################


function next{RKSF<:RKStepperFixed}(s::Solver{RKSF}, state::RKState)
    step = state.step
    work = state.work

    dof  = length(step.y)
    b    = s.stepper.tableau.b
    dt   = min(state.dt,s.options.tstop-step.t)

    # m3: why is it necessary to copy here and then copy back below?

    # pwl: to my understanding calc_next_k! needs the starting value
    # (step.y) but work.ynew is changing in the inner loop, so we need
    # two distinct arrays, both starting as step.y.

    copy!(work.ynew,step.y)

    for k=1:length(b)
        calc_next_k!(work, k, s.ode, s.stepper.tableau, step, dt)
        for d=1:dof
            work.ynew[d] += dt * b[k]*work.ks[k][d]
        end
    end
    step.t += dt
    copy!(step.y,work.ynew)
    return ((step.t,step.y), state)
end


########################
# Adaptive step method #
########################


function next{RKSA<:RKStepperAdaptive}(sol::Solver{RKSA}, state::RKState)

    const timeout_const = 5

    # the initial values
    dt      = state.dt          # dt is the previous stepisze, it is
    # modified inside the loop
    timeout = state.timeout
    work    = state.work
    step    = state.step
    tableau = sol.stepper.tableau

    # The while loop continues until we either find a stepsize which
    # leads to a small enough error or the stepsize reaches
    # prob.minstep

    # trim the inital stepsize to avoid overshooting
    dt      = min(dt, sol.options.tstop-state.step.t)

    while true

        # Do one step (assumes ks[1]==f0).  After calling work.ynew
        # holds the new step.
        # TODO: return ynew instead of passing it as work.ynew?

        # work.y and work.yerr and work.ks are updated after this step
        rk_embedded_step!(work, sol.ode, tableau, step, dt)

        # changes work.yerr
        err, newdt, timeout = stepsize_hw92!(work, step, tableau, dt, timeout, sol.options)

        # trim again in case newdt > dt
        newdt = min(newdt, sol.options.tstop-state.step.t)

        if abs(newdt) < sol.options.minstep  # minimum step size reached, break
            # passing the newdt to state will result in done()
            state.dt = newdt
            break
        end

        if err > 1 # error is too large, repeat the step with smaller dt
            # redo step with smaller dt and reset the timeout
            dt      = newdt
            timeout = timeout_const
        else
            # step is accepted

            # preload ks[1] for the next step
            if isFSAL(sol.stepper.tableau)
                copy!(work.ks[1],work.ks[end])
            else
                sol.ode.F!(step.t+dt, work.ynew, work.ks[1])
            end

            # Swap bindings of y and ytrial, avoids one copy
            step.y, work.ynew = work.ynew, step.y

            # Update state with the data from the step we have just
            # made:
            step.t += dt
            state.dt = newdt
            state.timeout = timeout
            break
        end
    end
    return ((step.t,step.y),state)
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

end


function stepsize_hw92!{T}(work,
                           last_step ::Step,
                           tableau   ::TableauRKExplicit,
                           dt        ::T,
                           timeout,
                           options   ::Options)
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
    fac = [T(8//10), T(9//10), T(1//4)^(1//(ord+1)), T(38//100)^(1//(ord+1))][1]
    facmax = T(5) # maximal step size increase. 1.5-5
    facmin = 1./facmax  # maximal step size decrease. ?
    dof = length(last_step.y)

    # in-place calculate yerr./tol
    for d=1:dof

        # if outside of domain (usually NaN) then make step size smaller by maximum
        if options.isoutofdomain(work.y[d])
            return T(10), dt*facmin, timout_after_nan
        end

        y0 = last_step.y[d] # TODO: is this supposed to be the last successful step?
        y1 = work.ynew[d]    # the approximation to the next step
        sci = (options.abstol + options.reltol*max(options.norm(y0),options.norm(y1)))
        work.yerr[d] = work.yerr[d]/sci # Eq 4.10
    end

    # TOOD: should we use options.norm here as well?
    err   = norm(work.yerr) # Eq. 4.11
    newdt = min(options.maxstep, dt*max(facmin, fac*(1/err)^(1//(ord+1)))) # Eq 4.13 modified

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
end
