function Step(problem :: RKProblem)
    t0 = problem.t0
    y0 = problem.y0
    dy0 = problem.F(t0,y0)
    dt0 = problem.dt0
    return Step(t0,y0,dy0,dt0)
end

type DenseState
    s0; s1
    last_tout
    first_step
    rkstate
    # used for storing the interpolation result
    ytmp
end

immutable DenseProblem
    rkprob :: RKProblem
    points :: Symbol
    tspan
end


function newDenseProblem(args...; tspan = [Inf], points = :all, opt_args...)
    rkprob = newRKProblem(args...; opt_args..., tstop = tspan[end])
    return DenseProblem(rkprob, points, tspan)
end

function start(prob :: DenseProblem)
    step0 = Step(prob.rkprob)
    step1 = Step(prob.rkprob)
    rkstate = start(prob.rkprob)
    ytmp = deepcopy(prob.rkprob.y0)
    return DenseState(step0, step1, prob.rkprob.t0, true, rkstate, ytmp)
end

function next(prob :: DenseProblem, state :: DenseState)

    s0, s1 = state.s0, state.s1
    t0, t1 = s0.t, s1.t

    if state.first_step
        state.first_step = false
        return ((s0.t,s0.y),state)
    end

    # the next output time that we aim at
    t_goal = prob.tspan[findfirst(t->(t>state.last_tout), prob.tspan)]

    # the t0 == t1 part ensures that we make at least one step
    while t1 < t_goal

        # s1 is the starting point for the new step, while the new
        # step is saved in s0

        if done(prob.rkprob, state.rkstate)
            error("The iterator was exhausted before the dense output compltede.")
        else
            # at this point s0 holds the new step, "s2" if you will
            ((s0.t,s0.y[:]),state.rkstate) = next(prob.rkprob, state.rkstate)
        end

        # swap s0 and s1
        s0, s1 = s1, s0
        # update the state
        state.s0, state.s1 = s0, s1
        # and times
        t0, t1 = s0.t, s1.t

        # we made a successfull step and points == :all
        if prob.points == :all
            t_goal = min(t_goal,t1)
            break
        end
    end

    # at this point we have t_goal∈[t0,t1] so we can apply the
    # interpolation

    F = prob.rkprob.F
    s0.dy[:], s1.dy[:] = F(t0,s0.y), F(t1,s1.y)

    hermite_interp!(state.ytmp,t_goal,s0,s1)

    # update the last output time
    state.last_tout = t_goal

    return ((t_goal,state.ytmp),state)

end

function done(prob :: DenseProblem, state :: DenseState)
    return done(prob.rkprob, state.rkstate) || state.s1.t >= prob.tspan[end]
end
