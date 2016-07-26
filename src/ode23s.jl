# ODE23S  Solve stiff systems based on a modified Rosenbrock triple
# (also used by MATLAB's ODE23s); see Sec. 4.1 in
#
# [SR97] L.F. Shampine and M.W. Reichelt: "The MATLAB ODE Suite," SIAM Journal on Scientific Computing, Vol. 18, 1997, pp. 1–22

immutable ModifiedRosenbrockIntegrator{T<:Number} <: AbstractIntegrator
    options::AdaptiveOptions{T}
    const_d::T
    const_e::T
end

@compat function (::Type{ModifiedRosenbrockIntegrator{T}}){T}(;options...)
    const_d = 1/(2+sqrt(T(2)))
    const_e = 6+sqrt(T(2))

    ModifiedRosenbrockIntegrator( AdaptiveOptions{T}(;options...), const_d, const_e )
end

order(::ModifiedRosenbrockIntegrator) = 2
name(::ModifiedRosenbrockIntegrator) = "Modified Rosenbrock Integrator"
isadaptive(::ModifiedRosenbrockIntegrator) = true
tdir(ode::ExplicitODE, stepper::ModifiedRosenbrockIntegrator) = sign(stepper.options.tstop - ode.t0)

# define the set of ODE problems with which this stepper can work
solve{T,S<:ModifiedRosenbrockIntegrator}(ode::ExplicitODE{T}, stepper::Type{S}; options...) =
    Problem(ode,stepper{T}(;options...))

"""
The state for the Rosenbrock stepper

- step:  Last successful step
- F1,F2: Work arrays for storing the intermediate values of y'
- J:     Jacobian
- iters: Number of successful steps made

"""
type RosenbrockState{T,Y} <: AbstractState
    step ::Step{T,Vector{Y}}
    dt   ::T
    F1   ::Vector{Y}
    F2   ::Vector{Y}
    k1   ::Vector{Y}
    k2   ::Vector{Y}
    k3   ::Vector{Y}
    ynew ::Vector{Y}
    dtold::T
    J    ::Matrix{Y}
    iters::Int
end


# for debugging
function show(io::IO, state::RosenbrockState)
    show(io,state.step)
    println("dt =$(state.dt)")
    println("F1 =$(state.F1)")
    println("F2 =$(state.F2)")
    println("J  =$(state.J)")
end


function init{T}(ode::ExplicitODE{T},
                 stepper::ModifiedRosenbrockIntegrator)
    t  = ode.t0
    dt = stepper.options.initstep
    y  = ode.y0
    dy = zero(y)

    J  = Array(eltype(y),length(y),length(y))

    step  = Step(t,copy(y),copy(dy))
    state = RosenbrockState(step,
                            dt,
                            zero(y), # F1
                            zero(y), # F2
                            zero(y), # k1
                            zero(y), # k2
                            zero(y), # k3
                            zero(y), # ynew
                            dt*0,    # dtnew
                            J,       # J
                            0)       # iters

    # initialize the derivative and the Jacobian
    ode.F!(t,y,step.dy)
    ode.J!(t,y,state.J)

    return state
end


function trialstep!(ode::ExplicitODE,
                    stepper::ModifiedRosenbrockIntegrator,
                    state::RosenbrockState)
    # unpack
    step    = state.step
    opts    = stepper.options
    F1, F2, J = state.F1, state.F2, state.J
    k1,k2,k3,ynew =  state.k1, state.k2, state.k3, state.ynew
    t, dt, y, dy = step.t, state.dt, step.y, step.dy
    F! = ode.F!
    F0 = dy

    td = tdir(ode,stepper)

    # see whether we're done
    if td*t >= td*opts.tstop
        # nothing left to integrate
        return finish
    end

    # increment iteration counter
    state.iters += 1
    if state.iters > opts.maxiters
        warn("Reached maximum number of iterations $(opts.maxiters)")
        return abort
    end

    W = lufact!( eye(J) - dt*stepper.const_d*J )

    # Approximate time-derivative of F, we are using F1 as a
    # temporary array
    F!(t+dt/100,y,F1)
    tder = 100*stepper.const_d*(F1-F0)

    # modified Rosenbrock formula
    # TODO: update k1,k2,k3 in-place
    k1[:] = W \ (F0 + tder)
    F!(t+dt/2, y+dt*k1/2, F1)
    k2[:] = W \ (F1 - k1) + k1
    for i=1:length(y)
        ynew[i] = y[i] + dt*k2[i]
    end
    F!(t+dt,   ynew,      F2)
    k3[:] = W \ (F2 - stepper.const_e*(k2 - F1) - 2*(k1 - F0) + tder )

    return cont
end

function errorcontrol!(ode::ExplicitODE,
                       stepper::ModifiedRosenbrockIntegrator,
                       state::RosenbrockState)

    step    = state.step
    opts    = stepper.options
    k1,k2,k3 =  state.k1, state.k2, state.k3
    k1,k2,k3,ynew =  state.k1, state.k2, state.k3, state.ynew
    t, dt, y, dy = step.t, state.dt, step.y, step.dy

    td = tdir(ode,stepper)

    # allowable error
    delta = max(opts.reltol*max(opts.norm(y), opts.norm(ynew)),opts.abstol)

    # error estimate
    err = (abs(dt)/6)*(opts.norm(k1 - 2*k2 + k3))/delta

    # new step-size
    dtnew = td*min(opts.maxstep, abs(dt)*0.8*err^(-1/3) )

    # trim in case newdt > dt
    dtnew = td*min(abs(dtnew), abs(opts.tstop-(t+dt)))

    state.dtold = dt
    state.dt = dtnew
    return err, cont
end

function accept!(ode::ExplicitODE,
                 stepper::ModifiedRosenbrockIntegrator,
                 state::RosenbrockState)
    step = state.step
    # update the state
    step.t     = step.t+state.dtold
    copy!(step.y, state.ynew)
    copy!(step.dy, state.F2)
    ode.J!(step.t,step.y,state.J)

    return cont
end
