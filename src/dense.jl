# A higher level stepper, defined as a wrapper around another stepper.

#TODO: how about having an DenseStepper <: AbstractWrapper <: AbstractStepper?
immutable DenseStepper <: AbstractStepper
    solver::Solver
end

solve(ode::ExplicitODE,
      stepper::DenseStepper,
      options::Options) = Solver(ode,stepper,options)

dense(sol::Solver) = solve(sol.ode, DenseStepper(sol), sol.options)

"""

The state of the dense stepper

- s0, s1: Previous steps, used to produce interpolated output
- solver_state: The state of the associated solver
- ytmp: work array

"""
type DenseState{St<:AbstractState,T,Y} <: AbstractState{T,Y}
    s0::Step{T,Y}
    s1::Step{T,Y}
    last_tout::T
    first_step
    solver_state::St
    # used for storing the interpolation result
    ytmp::Y
    solver_done
end


function start{O<:ExplicitODE,S<:DenseStepper}(s::Solver{O,S})
    # extract the real solver
    solver = s.stepper.solver
    t0  = solver.ode.t0
    y0  = solver.ode.y0
    dy0 = copy(y0)
    solver.ode.F!(t0,y0,dy0)
    step0 = Step(t0,copy(y0),copy(dy0))
    step1 = Step(t0,copy(y0),copy(dy0))
    solver_state = start(solver)
    ytmp = copy(y0)
    return DenseState(step0, step1, t0-1, true, solver_state, ytmp, false)
end

# m3: I think it would be nice to factor out the dense-output and
# root-finding into its own function.  That way it could be used also
# independently of the dense-output iterator.  Also, it would make
# this next function more compact.

# pwl: I agree, but then the problem is that once you decouple them
# you would lose the opprotunity to detect the roots with each step.
function next{O<:ExplicitODE,S<:DenseStepper}(s::Solver{O,S}, state::DenseState)


    # m3: I'm not 100% sure what happens here.  I would implement it like so:
    # Initialize in `start`: calculate next t1, y1 and also hold onto IC in t0,y0,
    # set state.last_t=1
    #
    # in next have this loop
    # for t in tspan[state.last_t:end]
    #     if t>t1
    #        make new t1, y1, move old t1, y1 into t0, y0
    #     end
    #     make dense output at t
    #     find events
    #     state.last_t += 1
    #     return ((t, y), state)
    # end

    # pwl: @m3 this is basically what happens here:-), although I'm
    # not using the index of tspan anywhere explicitly.

    solver = s.stepper.solver

    # these guys store the intermediate steps we make
    s0, s1 = state.s0, state.s1
    t0, t1 = s0.t, s1.t

    # assuming the last output was done at state.last_tout set the
    # t_goal to the next larger time from tspan.  Strong inequality
    # below is crucial, otherwise we would be selecting the same step
    # every time.
    tspan  = s.options.tspan
    t_goal = tspan[findfirst(t->(t>state.last_tout), tspan)]

    # Keep computing new steps (i.e. new pairs (t0,t1)) until we reach
    # t0 < t_goal <= t1, then we use interpolation to get the value at
    # t_goal.  Unless points==:all, then we break the while loop after
    # making the first step.
    while t_goal > t1

        # s1 stores the last succesfull step, the new step is stored
        # in s0

        if done(solver, state.solver_state)
            warn("The iterator was exhausted before the dense output completed.")
            # prevents calling done(..) twice
            state.solver_done = true
            # TODO: deepcopy?
            # Return whatever we got as the last step
            return ((s0.t,s0.y[:]),state)
        else
            # at this point s0 is updated with the new step, "s2" if you will
            ((s0.t,s0.y[:]), state.solver_state) = next(solver, state.solver_state)
        end

        # swap s0 and s1
        s0, s1 = s1, s0
        # and times
        t0, t1 = s0.t, s1.t

        # update the state accordingly
        state.s0, state.s1 = s0, s1

        # we haven't reached t_goal yet (t1<t_goal) but the option
        # points==:all calls for an output at every successful
        # step.
        if s.options.points == :all && t1 < t_goal-eps(t1)
            state.last_tout = t1
            return ((t1,s1.y),state)
        end

        if s.options.stopevent(t1,s1.y)
            break
        end

    end

    # at this point we have t0 < t_goal < t1 so we can apply the
    # interpolation to get a value of the solution at t_goal

    # TODO: is this necessary?  The solver should store the value of dy.
    solver.ode.F!(t0,s0.y,s0.dy)
    solver.ode.F!(t1,s1.y,s1.dy)

    if s.options.stopevent(t1,s1.y)
        function stopfun(t)
            hermite_interp!(state.ytmp,t,s0,s1)
            res = Int(s.options.stopevent(t,state.ytmp))
            return 2*res-1      # -1 if false, +1 if true
        end
        t_goal = findroot(stopfun, [s0.t,s1.t], s.options.roottol)
        # state.ytmp is already overwriten to the correct result as a
        # side-effect of calling stopfun
    else
        hermite_interp!(state.ytmp,t_goal,s0,s1)
    end

    # update the last output time
    state.last_tout = t_goal

    return ((t_goal,state.ytmp),state)

end


function done{O<:ExplicitODE,S<:DenseStepper}(s::Solver{O,S}, state::DenseState)

    return (
            state.solver_done ||
            state.last_tout >= s.options.tspan[end] ||
            s.options.stopevent(state.s1.t,state.s1.y)
            )
end


function hermite_interp!(y,t,step0::Step,step1::Step)
    # For dense output see Hairer & Wanner p.190 using Hermite
    # interpolation. Updates y in-place.
    #
    # f_0 = f(x_0 , y_0) , f_1 = f(x_0 + h, y_1 )
    # this is O(3). TODO for higher order.

    y0,  y1  = step0.y, step1.y
    dy0, dy1 = step0.dy, step1.dy

    if t == step0.t
        copy!(y,y0)
    elseif t == step1.t
        copy!(y,y1)
    else
        dt       = step1.t-step0.t
        theta    = (t-step0.t)/dt
        for i=1:length(y0)
            y[i] = ((1-theta)*y0[i] + theta*y1[i] + theta*(theta-1) *
                    ((1-2*theta)*(y1[i]-y0[i]) + (theta-1)*dt*dy0[i] + theta*dt*dy1[i]) )
        end
    end
end
