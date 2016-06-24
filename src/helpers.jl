#TODO make it a function on ExplicitODE and Options

"""

Chooses an initial step-size basing on the equation, initial data,
time span and the order of the method of integration.

"""
function dtinit{T}(F, y0, tspan::Vector{T}, reltol, abstol; order = 1)
    t0 = abs(tspan[1])
    tstop = abs(tspan[end])
    tau = max(reltol*norm(y0, Inf), abstol)
    d0 = norm(y0, Inf)/tau
    f0 = F(t0, y0)
    d1 = norm(f0, Inf)/tau
    if min(d0,d1) < eps(T)^(1/3)
        dt0 = eps(T)^(1/3)/10
    else
        dt0 = (d0/d1)/100
    end
    # perform Euler step
    y1 = y0+dt0*f0
    f1 = F(t0 + dt0, y1)
    # estimate second derivative
    d2 = norm(f1 - f0, Inf)/(tau*dt0)
    if max(d1, d2) <= 10*eps(T)
        dt1 = max(eps(T)^(1/3)/10, dt0/10^3)
    else
        pow = -(2 + log10(max(d1, d2)))/(order+1)
        dt1 = 10^pow
    end
    return min(100*dt0, dt1, abs(tstop-t0))
end

# a scalar version of the above
dtinit(F, y0::Number, args...; kargs...) = dtinit((t,y)->[F(t,y[1])], [y0], args...; kargs...)

"""

A simple bisection algorithm for finding a root of a solution f(x)=0
starting within the range x∈rng, the result is a point x₀ which is
located within the distance eps from the true root of f(x)=0.  For
this algorithm to work we need f(rng[1]) to have a different sign then
f(rng[2]).

"""
function findroot(f,rng,eps)
    xl, xr = rng
    fl, fr = f(xl), f(xr)

    if fl*fr > 0 || xl > xr
        error("Inconsistent bracket")
    end

    while xr-xl > eps
        xm = (xl+xr)/2
        fm = f(xm)

        if fm*fr > 0
            xr = xm
            fr = fm
        else
            xl = xm
            fl = fm
        end
    end

    return (xr+xl)/2
end
