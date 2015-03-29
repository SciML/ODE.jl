#ODE23  Solve non-stiff differential equations.
#
#   ODE23(F,TSPAN,Y0) with TSPAN = [T0 TFINAL] integrates the system
#   of differential equations dy/dt = f(t,y) from t = T0 to t = TFINAL.
#   The initial condition is y(T0) = Y0.
#
#   The first argument, F, is a function handle or an anonymous function
#   that defines f(t,y).  This function must have two input arguments,
#   t and y, and must return a column vector of the derivatives, dy/dt.
#
#   With two output arguments, [T,Y] = ODE23(...) returns a column
#   vector T and an array Y where Y(:,k) is the solution at T(k).
#
#   More than four input arguments, ODE23(F,TSPAN,Y0,RTOL,P1,P2,...),
#   are passed on to F, F(T,Y,P1,P2,...).
#
#   ODE23 uses the Runge-Kutta (2,3) method of Bogacki and Shampine (BS23).
#
#   Example
#      tspan = [0, 2*pi]
#      y0 = [1, 0]
#      F = (t, y) -> [0 1; -1 0]*y
#      ode23(F, tspan, y0)
#
# Adapted from Cleve Moler's textbook
# http://www.mathworks.com/moler/ncm/ode23tx.m
function ode23_ncm(F, y0, tspan; reltol = 1.e-5, abstol = 1.e-8)
    if reltol == 0
        warn("setting reltol = 0 gives a step size of zero")
    end

    threshold = abstol / reltol

    t = tspan[1]
    tfinal = tspan[end]
    tdir = sign(tfinal - t)
    hmax = abs(0.1*(tfinal-t))
    y = y0

    tout = Array(typeof(t), 1)
    tout[1] = t         # first output time
    yout = Array(typeof(y0),1)
    yout[1] = y         # first output solution

    # Compute initial step size.
    s1 = F(t, y)
    r = norm(s1./max(abs(y), threshold), Inf) + realmin() # TODO: fix type bug in max()
    h = tdir*0.8*reltol^(1/3)/r

    # The main loop.

    while t != tfinal

        hmin = 16*eps()*abs(t)
        if abs(h) > hmax; h = tdir*hmax; end
        if abs(h) < hmin; h = tdir*hmin; end

        # Stretch the step if t is close to tfinal.

        if 1.1*abs(h) >= abs(tfinal - t)
            h = tfinal - t
        end

        # Attempt a step.

        s2 = F(t+h/2, y+h/2*s1)
        s3 = F(t+3*h/4, y+3*h/4*s2)
        tnew = t + h
        ynew = y + h*(2*s1 + 3*s2 + 4*s3)/9
        s4 = F(tnew, ynew)

        # Estimate the error.

        e = h*(-5*s1 + 6*s2 + 8*s3 - 9*s4)/72
        err = norm(e./max(max(abs(y), abs(ynew)), threshold), Inf) + realmin()

        # Accept the solution if the estimated error is less than the tolerance.

        if err <= reltol
            t = tnew
            y = ynew
            push!(tout, t)
            push!(yout, y)
            s1 = s4   # Reuse final function value to start new step
        end

        # Compute a new step size.

        h = h*min(5, 0.8*(reltol/err)^(1/3))

        # Exit early if step size is too small.

        if abs(h) <= hmin
            println("Step size ", h, " too small at t = ", t)
            t = tfinal
        end

    end # while (t != tfinal)

    return tout, yout

end # ode23_ncm
