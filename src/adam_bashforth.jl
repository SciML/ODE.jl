####################################
# Explicit Adam-Bashforth solvers
####################################
# (Hairer & Wanner 1996, Vol I, p.357-358
#Begin AdamBash module
module AdamBash

using Polynomials;
using Compat;
include("ODE.jl")

##Implementation of the Adam Steps for order 2, 3 and 4 using recursive formualae from Hairer
abstract Integrator
immutable AdamBashforth{N} end

@inline function adam_step(::Union{AdamBashforth{2},Type{AdamBashforth{2}}}, ydot, y, stepsize, i)
    ynext = y[i] + stepsize/2*(3*ydot[i]
        - 1*ydot[i-1])
end
@inline function adam_step(::Union{AdamBashforth{3},Type{AdamBashforth{3}}}, ydot, y, stepsize, i)
    ynext = y[i] + stepsize/12*(23*ydot[i]
    - 16*ydot[i-1]
    + 5*ydot[i-2])
end
@inline function adam_step(::Union{AdamBashforth{4},Type{AdamBashforth{4}}}, ydot, y, stepsize, i)
    ynext = y[i] + stepsize/24*(55*ydot[i]
    - 59*ydot[i-1]
    + 37*ydot[i-2]
    - 9*ydot[i-3])
end

##Function: ode_abe(F,y0, tspan, order(optional))
##order is set to 4 by default
##Adam Bashforth is a fixed step, fixed order, multistep explicit method
function ode_ab(F::Function,y0, tspan,order=4 ::Integer)
    if !(2 <= order <= 4)
        error("Currently only orders 2,3, and 4 are implemented for Adam-Bashforth method")
    end

    h = diff(tspan)
    y = Array(typeof(y0), length(tspan))
    ydot = similar(y)
    y[1] = y0
    ydot[1] = F(tspan[1],y[1])

    ##Use Runge-Kunta for initial points necessary to base Adam Basforth off of, if initial values not given
    tint, yint= ODE.ode_ms(F,y0, tspan[1:order], order)
    for i = 1 : order
        y[i] = yint[i]
        ydot[i] = F(tspan[i],y[i])
    end
    ##Use Adam Basforth method for subsequent steps
    Integrator = AdamBashforth{order}
    for i = order+1:length(tspan)
        y[i] = adam_step(Integrator,ydot,y,h[i-1],i-1)
        ydot[i] = F(tspan[i],y[i])
    end
    ##Return y values
    return tspan, y
end
end #End Module
