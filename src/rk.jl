# This file contains the implementation of explicit Runkge-Kutta
# solver from (Hairer & Wanner 1992 p.134, p.165-169).

using Iterators

# include the Butcher tableaus.
include("tableaus.jl")

####################
# Iterator methods #
####################

# common structures and functions

type TempArrays
    y; ks; yerr
end

type State
    t; dt; y; tmp :: TempArrays
    timeout
end


immutable Problem{MethodType}
    method
    F
    y0
    t0
    dt0
    tstop
    reltol
    abstol
    minstep
    maxstep
end


# overload the call method for TableauRKExplicit, it returns the
# iterator (fixed or variable step according to the tableau)
function call(tab::TableauRKExplicit,
              F, y0, t0;
              tstop  = Inf,
              reltol = 1e-5,
              abstol = 1e-5,
              minstep = 1e-10,
              maxstep = 1/minstep,
              dt0 = hinit(F, y0, t0, tstop, tab, reltol, abstol),
              kargs...
              )

    if isadaptive(tab)
        methodtype = :adaptive
    else
        methodtype = :fixed
    end

    return Problem{methodtype}(tab, F, y0, t0, dt0, tstop, reltol, abstol, minstep, maxstep)
end


function start(problem :: Problem)
    t0, dt0, y0 = problem.t0, problem.dt0, problem.y0

    tmp = TempArrays(similar(y0), Array(typeof(y0), S(problem.method)), similar(y0))
    tmp.ks[1] = problem.F(t0,y0) # we assume that ks[1] is already initialized

    timeout = 0 # for step control
    return State(t0,dt0,deepcopy(y0),tmp,timeout)
end


function done(prob :: Problem, state :: State)
    return state.t >= prob.tstop || state.dt < prob.minstep
end


#####################
# Fixed step method #
#####################


function next(prob :: Problem{:fixed}, state :: State)
    dof = length(state.y)
    for s=1:S(prob.method)
        calc_next_k!(state.tmp, s, state, prob)
        for d=1:dof
            state.y[d] += state.dt * prob.method.b[s]*state.tmp.ks[s][d]
        end
    end
    state.t += state.dt
    return ((state.t,state.y), state)
end


########################
# Adaptive step method #
########################


function next(prob :: Problem{:adaptive}, state :: State)

    const timeout_const = 5

    # the initial values
    dt      = state.dt          # dt is the previous stepisze, it is
                                # modified inside the loop
    timeout = state.timeout

    # for aesthetical reasons we extract the temporary componen
    tmp     = state.tmp

    # The while loop continues until we either find a stepsize which
    # leads to a small enough error or the stepsize reaches
    # prob.minstep

    while true

        # do one step (assumes ks[1]==f0), changes only tmp
        err, newdt, timeout = rk_trial_step!(tmp, state, prob, dt, timeout)

        if abs(newdt) < prob.minstep  # minimum step size reached, break
            println("Warning: dt < minstep.  Stopping.")
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
            if isFSAL(prob.method)
                tmp.ks[1] = tmp.ks[S(prob.method)]
            else
                tmp.ks[1] = prob.F(state.t+dt, state.tmp.y)
            end

            # Swap bindings of y and ytrial, avoids one copy
            state.y, state.tmp.y = state.tmp.y, state.y

            # Update state with the data from the step we have just
            # made:
            state.t += dt
            state.dt = newdt
            state.timeout = timeout
            break
        end
    end
    return ((state.t,state.y),state)
end


##########################
# Lower level algorithms #
##########################


function rk_trial_step!(tmp, state, prob, dt, timeout)

    # tmp.y and tmp.yerr and tmp.ks are updated after this step
    rk_embedded_step!(tmp, state, prob, dt)

    # changes tmp.yerr (via in place update)
    err, newdt, timeout = stepsize_hw92!(tmp, state, prob, dt, timeout)

    return err, newdt, timeout
end


function rk_embedded_step!(tmp :: TempArrays, state :: State, prob :: Problem, dt)
    # Does one embedded R-K step updating ytrial, yerr and ks.
    #
    # Assumes that ks[:,1] is already calculated!
    #
    # Modifies tmp.y and tmp.yerr only

    y      = state.y
    dof    = length(y)
    b      = prob.method.b

    tmp.y[:]    = 0
    tmp.yerr[:] = 0

    for d=1:dof
        tmp.y[d]    += b[1,1]*tmp.ks[1][d]
        tmp.yerr[d] += b[2,1]*tmp.ks[1][d]
    end

    for s=2:S(prob.method)
        calc_next_k!(state.tmp, s, state, prob)
        for d=1:dof
            tmp.y[d]    += b[1,s]*tmp.ks[s][d]
            tmp.yerr[d] += b[2,s]*tmp.ks[s][d]
        end
    end

    for d=1:dof
        tmp.yerr[d] = dt * (tmp.y[d]-tmp.yerr[d])
        tmp.y[d]    = y[d] + dt * tmp.y[d]
    end
end


function stepsize_hw92!(tmp, state, prob, dt, timeout)
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

    order = minimum(prob.method.order)
    timout_after_nan = 5
    fac = [0.8, 0.9, 0.25^(1/(order+1)), 0.38^(1/(order+1))][1]
    facmax = 5.0 # maximal step size increase. 1.5-5
    facmin = 1./facmax  # maximal step size decrease. ?
    dof = length(state.y)

    # in-place calculate yerr./tol
    for d=1:dof

        # if outside of domain (usually NaN) then make step size smaller by maximum
        if isoutofdomain(tmp.y[d])
            return 10., dt*facmin, timout_after_nan
        end

        tmp.yerr[d] = tmp.yerr[d]/(prob.abstol + max(norm(prob.y0[d]), norm(tmp.y[d]))*prob.reltol) # Eq 4.10
    end

    err = norm(tmp.yerr, 2) # Eq. 4.11
    newdt = min(prob.maxstep, dt*max(facmin, fac*(1/err)^(1/(order+1)))) # Eq 4.13 modified

    if timeout > 0
        newdt = min(newdt, dt)
        timeout -= 1
    end

    return err, newdt, timeout
end


function hinit(F, y0, t0, tstop, method, reltol, abstol)
    # Returns first step size
    tdir = sign(tstop - t0)
    order = minimum(method.order)
    tau = max(reltol*norm(y0, Inf), abstol)
    d0 = norm(y0, Inf)/tau
    f0 = F(t0, y0)
    d1 = norm(f0, Inf)/tau
    if d0 < 1e-5 || d1 < 1e-5
        h0 = 1e-6
    else
        h0 = 0.01*(d0/d1)
    end
    # perform Euler step
    y1 = y0 + tdir*h0*f0
    f1 = F(t0 + tdir*h0, y1)
    # estimate second derivative
    d2 = norm(f1 - f0, Inf)/(tau*h0)
    if max(d1, d2) <= 1e-15
        h1 = max(1e-6, 1e-3*h0)
    else
        pow = -(2 + log10(max(d1, d2)))/(order+1)
        h1 = 10^pow
    end
    return min(100*h0, h1, tdir*abs(tstop-t0))
end


# For clarity we pass the TempArrays part of the state separately,
# this is the only part of state that can be changed here
function calc_next_k!(tmp :: TempArrays, i, state :: State, prob :: Problem)
    dof = length(state.y)
    t, dt, a, c = state.t, state.dt, prob.method.a, prob.method.c

    tmp.y[:] = state.y
    for j=1:i-1
        for d=1:dof
            tmp.y[d] += dt * tmp.ks[j][d] * a[i,j]
        end
    end
    tmp.ks[i] = prob.F(t + c[i]*dt, tmp.y)

    nothing
end
