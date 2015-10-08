####################
# Iterator methods #
####################

type TempArrays
    y; ks
end

type RKState
    t; dt; y; tmp :: TempArrays
end

immutable RKProblem
    F
    btab
    y0
    t0
    dt0
    tstop
end

function newRKProblem(fn, y0, t0, dt0; tstop = Inf, method = bt_feuler)
    return RKProblem(fn, method, y0, t0, dt0, tstop)
end

function start(problem :: RKProblem)
    t0  = problem.t0
    dt0 = problem.dt0
    y0  = problem.y0
    tmp = TempArrays(problem.y0, Array(typeof(y0), S(problem.btab)))
    return RKState(t0, dt0, y0, tmp)
end

function next(prob :: RKProblem, state :: RKState)

    dof = length(state.y)
    for s=1:S(prob.btab)
        calc_next_k!(state.tmp, s, state, prob)
        for d=1:dof
            state.y[d] += state.dt * prob.btab.b[s]*state.tmp.ks[s][d]
        end
    end

    state.t += state.dt

    return ((state.t,state.y), state)

end

function done(prob :: RKProblem, state :: RKState)
    return state.t >= prob.tstop
end


function calc_next_k!(tmp :: TempArrays, i, state :: RKState, prob :: RKProblem)
    dof = length(state.y)
    t, dt, a, c = state.t, state.dt, prob.btab.a, prob.btab.c

    tmp.y[:] = state.y
    for j=1:i-1
        # tmp.y += dt * btab.a[i,j] *  ks[j]
        for d=1:dof
            tmp.y[d] += dt * tmp.ks[j][d] * a[i,j]
        end
    end
    tmp.ks[i] = prob.F(t + c[i]*dt, tmp.y)

    nothing
end
