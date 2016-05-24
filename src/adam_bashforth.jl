####################################
# Explicit Adam-Bashforth solvers
####################################
# (Hairer & Wanner 1996, Vol I, p.357-358
using Polynomials;
using Compat;
include("ODE.jl")

##Implementation of the Adam Steps for order 2, 3 and 4
abstract Integrator
immutable AdamBashforth{N} end

@inline function adam_step(::Union{AdamBashforth{2},Type{AdamBashforth{2}}}, f, y, t, stepsize, i)
    ynext = y[i] + stepsize/2*(3*f(t[i], y[i])
    - 1*f(t[i-1],y[i-1]))
end
@inline function adam_step(::Union{AdamBashforth{3},Type{AdamBashforth{3}}}, f, y, t, stepsize, i)
    ynext = y[i] + stepsize/12*(23*f(t[i], y[i])
    - 16*f(t[i-1],y[i-1])
    + 5*f(t[i-2],y[i-2]))
end
@inline function adam_step(::Union{AdamBashforth{4},Type{AdamBashforth{4}}}, f, y, t, stepsize, i)
    ynext = y[i] + stepsize/24*(55*f(t[i], y[i])
    - 59*f(t[i-1],y[i-1])
    + 37*f(t[i-2],y[i-2])
    - 9*f(t[i-3],y[i-3]))
end

##Fixed step, fixed order, multistep explicit method
function ode_ab(F,y0, tspan,order ::Integer)
    if !(2 <= order <= 4)
        error("Currently only orders 2,3, and 4 are implemented for Adam-Bashforth method")
    end

    h = diff(tspan)
    yvals = Array(typeof(y0), length(tspan))
    yvals[1] = y0

    ##Use Runge-Kunta for initial points necessary to base Adam Basforth off of, if initial values not given
    tint, yint= ODE.ode_ms(F,y0, tspan[1:order], order)

    ##Use Adam Basforth method for subsequent steps
    for i = 1 : order
        yvals[i] = yint[i]
    end
    Integrator = AdamBashforth{order}
    for i = order+1:length(tspan)
        yvals[i] = adam_step(Integrator,F,yvals,tspan,h[i-1],i-1)
        #yvals[i] = y
    end
    ##Return y values
    tspan, yvals
end
