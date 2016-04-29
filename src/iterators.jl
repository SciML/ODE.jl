"""
Generic done method, some steppers may implement their own versions.
"""


function done(s::Solver, state::AbstractState)
    if state.step.t >= s.options.tstop
        return true
    elseif state.dt < s.options.minstep
        warn("minstep reached.")
        return true
    elseif state.iters >= s.options.maxiters
        warn("Maximum number of iterations ($(Int(s.options.maxiters))) reached, consider setting a larger maxiter.")
        return true
    end
    return false
end
