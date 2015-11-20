import Base: start, next, done, call, collect

include("dense.jl")
include("rk.jl")

# this wraps any iterator (method) returning pairs (t,y) in a dense
# output and also covers the reverse time integration
function solver(F, y0::AbstractArray, t0;
                tstop = Inf,
                tspan = [tstop],
                method = bt_feuler,
                stopevent = (t,y)->false,
                kargs...)

    if tstop >= t0
        # forward time integration
        sol = method(F,y0,t0; tstop = tstop, kargs...)
        dense_sol = dense(F, y0, t0, sol; tspan = tspan, kargs...)
        return dense_sol
    else
        # reverse time integration
        F_reverse(t,y) = -F(2*t0-t,y)
        reverse_output(t,y)=(2*t0-t,y)
        sol = method(F_reverse,y0,t0;
                     tstop = 2*t0 -tstop,
                     tspan = 2*t0.-tspan,
                     kargs...)
        dense_sol = dense(F_reverse, y0, t0, sol;
                          tspan = 2*t0-tspan,
                          stopevent = (t,y)->stopevent(2*t0-t,y),
                          kargs...)

        return imap(x->reverse_output(x...),dense_sol)
    end

end


solver(F,y0,t0;kargs...)=solver((t,y)->[F(t,y[1])],[y0],t0;kargs...)
