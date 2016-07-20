"""
Generic done method, some steppers may implement their own versions.
"""

function done(s::Solver, state::AbstractState)
    st = s.stepper
    if state.step.t >= st.options.tstop
        return true
    elseif state.dt < st.options.minstep
        warn("minstep reached.")
        return true
    elseif state.iters >= st.options.maxiters
        warn("Maximum number of iterations ($(Int(s.options.maxiters))) reached, consider setting a larger maxiter.")
        return true
    end
    return false
end
