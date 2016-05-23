####################################
# Explicit Adam-Bashforth solvers
####################################
# (Hairer & Wanner 1996, Vol I, p.357-358

include("runge_kutta.jl")

##Fixed step, fixed order, multistep explicit method 

function odeadambashforth{N}(F,y0, tspan,order ::Integer)
    if !(2 <= order <= 4)
      error("Currently only orders 2,3, and 4 are implemented for Adam-Bashforth method")
    
    h = diff(tspan)
    y = Array(typeof(y0), length(tspan))
    y[1] = x0
    
    ##Use Runge-Kunta for initial points necessary to base Adam Basforth off of
    @inbounds for i = 0:N-1
        y[i] = use runge-kunta
    end
    
    ##Use Adam Basforth method for subsequent steps
    Integrator = AdamBashforth{order}
    @inbounds for i = N:niter
        y = step(Integrator,f,y,t,stepsize,i-1)
        y[i] = y
        t = t + stepsize
    end
    ##Return y values
    yvals
end

##Implementation of the Adam Steps for order 2, 3 and 4

abstract Integrator
immutable AdamBashforthStep{N} end

@inline function adam_step(::Union{AdamBashforthStep{2},Type{AdamBashforthStep{2}}}, f, y, t, stepsize)
    ynext = y[i] + stepsize/2*(3*f(t[i], y[i]) 
                                - 1*f(t[i-1],y[i-1]))
end

@inline function step(::Union{AdamBashforthStep{3},Type{AdamBashforthStep{3}}}, f, y, t, stepsize)
      ynext = y[i] + stepsize/12*(23*f(t[i], y[i]) 
                                - 16*f(t[i-1],y[i-1])
                                + 5*f(t[i-2],y[i-2]))
end


@inline function step(::Union{AdamBashforthStep{4},Type{AdamBashforthStep{4}}}, f, y, t, stepsize)
      ynext = y[i] + stepsize/24*(55*f(t[i], y[i]) 
                                - 59*f(t[i-1],y[i-1])
                                + 37*f(t[i-2],y[i-2])
                                - 9*f(t[i-3],y[i-3]))
end
