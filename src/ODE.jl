# Ordinary Differential Equation Solvers

module ODE

using Polynomials

## minimal function export list
# adaptive non-stiff:
export ode23, ode45, ode78
# adaptive stiff:
export ode23s
# non-adaptive:
export ode4s, ode4ms, ode4

## complete function export list
#export ode23, ode4,
#    oderkf, ode45, ode45_dp, ode45_fb, ode45_ck,
#    oderosenbrock, ode4s, ode4s_kr, ode4s_s,
#    ode4ms, ode_ms

###############################################################################
## HELPER FUNCTIONS
###############################################################################

# estimator for initial step based on book
# "Solving Ordinary Differential Equations I" by Hairer et al., p.169
function hinit(F, x0, t0, tdir, p, reltol, abstol)
    tau = max(reltol*norm(x0, Inf), abstol)
    d0 = norm(x0, Inf)/tau
    f0 = F(t0, x0)
    d1 = norm(f0, Inf)/tau
    if d0 < 1e-5 || d1 < 1e-5
        h0 = 1e-6
    else
        h0 = 0.01*(d0/d1)
    end
    # perform Euler step
    x1 = x0 + tdir*h0*f0
    f1 = F(t0 + tdir*h0, x1)
    # estimate second derivative
    d2 = norm(f1 - f0, Inf)/(tau*h0)
    if max(d1, d2) <= 1e-15
        h1 = max(1e-6, 1e-3*h0)
    else
        pow = -(2. + log10(max(d1, d2)))/(p + 1.)
        h1 = 10.^pow
    end
    h = min(100*h0, h1), f0
end

###############################################################################
## NON-STIFF SOLVERS
###############################################################################

# ODERKF based on
#
# ode45 adapted from http://users.powernet.co.uk/kienzle/octave/matcompat/scripts/ode_v1.11/ode45.m
# (a newer version (v1.15) can be found here https://sites.google.com/site/comperem/home/ode_solvers)
#
# Original Octave implementation:
# Marc Compere
# CompereM@asme.org
# created : 06 October 1999
# modified: 17 January 2001
#
function oderkf(F, x0, tspan, p, a, bs, bp; reltol = 1.0e-5, abstol = 1.0e-8,
                                            norm=Base.norm,
                                            minstep=abs(tspan[end] - tspan[1])/1e9,
                                            maxstep=abs(tspan[end] - tspan[1])/2.5,
                                            initstep=0.)
    # see p.91 in the Ascher & Petzold reference for more infomation.
    pow = 1/p   # use the higher order to estimate the next step size

    c = sum(a, 2)   # consistency condition
    k = Array(typeof(x0), length(c))

    # Initialization
    t = tspan[1]
    tfinal = tspan[end]
    tdir = sign(tfinal - t)

    h = initstep
    if h == 0.
      # initial guess at a step size
      h, k[1] = hinit(F, x0, t, tdir, p, reltol, abstol)
    else
      k[1] = F(t, x0) # first stage
    end
    h = tdir*min(h, maxstep)
    x = x0
    tout = Array(typeof(t), 1)
    tout[1] = t         # first output time
    xout = Array(typeof(x0), 1)
    xout[1] = x         # first output solution

    while abs(t) != abs(tfinal) && abs(h) >= minstep
        if abs(h) > abs(tfinal-t)
            h = tfinal - t
        end

        #(p-1)th and pth order estimates
        xs = x + h*bs[1]*k[1]
        xp = x + h*bp[1]*k[1]
        for j = 2:length(c)
            dx = a[j,1]*k[1]
            for i = 2:j-1
                dx += a[j,i]*k[i]
            end
            k[j] = F(t + h*c[j], x + h*dx)

            # compute the (p-1)th order estimate
            xs = xs + h*bs[j]*k[j]
            # compute the pth order estimate
            xp = xp + h*bp[j]*k[j]
        end

        # estimate the local truncation error
        gamma1 = xs - xp

        # Estimate the error and the acceptable error
        delta = norm(gamma1, Inf)              # actual error
        tau   = max(reltol*norm(x,Inf),abstol) # allowable error

        # Update the solution only if the error is acceptable
        if delta <= tau
            t = t + h
            x = xp    # <-- using the higher order estimate is called 'local extrapolation'
            push!(tout, t)
            push!(xout, x)

            # Compute the slopes by computing the k[:,j+1]'th column based on the previous k[:,1:j] columns
            # notes: k needs to end up as an Nxs, a is 7x6, which is s by (s-1),
            #        s is the number of intermediate RK stages on [t (t+h)] (Dormand-Prince has s=7 stages)
            if c[end] == 1
                # Assign the last stage for x(k) as the first stage for computing x[k+1].
                # This is part of the Dormand-Prince pair caveat.
                # k[:,7] has already been computed, so use it instead of recomputing it
                # again as k[:,1] during the next step.
                k[1] = k[end]
            else
                k[1] = F(t,x) # first stage
            end
        end

        # Update the step size
        h = tdir*min(maxstep, 0.8*abs(h)*(tau/delta)^pow)
    end # while (t < tfinal) & (h >= minstep)

    if abs(t) < abs(tfinal)
        error("Step size grew too small. t=$t, h=$(abs(h)), x=$x")
    end

    return tout, xout
end


# Bogacki–Shampine coefficients
const bs_coefficients = (3,
                         [    0           0      0      0
                              1/2         0      0      0
                              0         3/4      0      0
                              2/9       1/3     4/9     0],
                         # 2nd order b-coefficients
                         [7/24 1/4 1/3 1/8],
                         # 3rd order b-coefficients
                         [2/9 1/3 4/9 0],
                         )
ode23_bs(F, x0, tspan; kwargs...) = oderkf(F, x0, tspan, bs_coefficients...; kwargs...)


# Both the Dormand-Prince and Fehlberg 4(5) coefficients are from a tableau in
# U.M. Ascher, L.R. Petzold, Computer Methods for  Ordinary Differential Equations
# and Differential-Agebraic Equations, Society for Industrial and Applied Mathematics
# (SIAM), Philadelphia, 1998
#
# Dormand-Prince coefficients
const dp_coefficients = (5,
                         [    0           0          0         0         0        0
                              1/5         0          0         0         0        0
                              3/40        9/40       0         0         0        0
                             44/45      -56/15      32/9       0         0        0
                          19372/6561 -25360/2187 64448/6561 -212/729     0        0
                           9017/3168   -355/33   46732/5247   49/176 -5103/18656  0
                             35/384       0        500/1113  125/192 -2187/6784  11/84],
                         # 4th order b-coefficients
                         [5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40],
                         # 5th order b-coefficients
                         [35/384 0 500/1113 125/192 -2187/6784 11/84 0],
                         )
ode45_dp(F, x0, tspan; kwargs...) = oderkf(F, x0, tspan, dp_coefficients...; kwargs...)


# Fehlberg coefficients
const fb_coefficients = (5,
                         [    0         0          0         0        0
                             1/4        0          0         0        0
                             3/32       9/32       0         0        0
                          1932/2197 -7200/2197  7296/2197    0        0
                           439/216     -8       3680/513  -845/4104   0
                            -8/27       2      -3544/2565 1859/4104 -11/40],
                         # 4th order b-coefficients
                         [25/216 0 1408/2565 2197/4104 -1/5 0],
                         # 5th order b-coefficients
                         [16/135 0 6656/12825 28561/56430 -9/50 2/55],
                         )
ode45_fb(F, x0, tspan; kwargs...) = oderkf(F, x0, tspan, fb_coefficients...; kwargs...)


# Cash-Karp coefficients
# Numerical Recipes in Fortran 77
const ck_coefficients = (5,
                         [   0         0       0           0          0
                             1/5       0       0           0          0
                             3/40      9/40    0           0          0
                             3/10     -9/10    6/5         0          0
                           -11/54      5/2   -70/27       35/27       0
                          1631/55296 175/512 575/13824 44275/110592 253/4096],
                         # 4th order b-coefficients
                         [37/378 0 250/621 125/594 0 512/1771],
                         # 5th order b-coefficients
                         [2825/27648 0 18575/48384 13525/55296 277/14336 1/4],
                         )
ode45_ck(F, x0, tspan; kwargs...) = oderkf(F, x0, tspan, ck_coefficients...; kwargs...)


# Fehlberg 7(8) coefficients
# Values from pag. 65, Fehlberg, Erwin. "Classical fifth-, sixth-, seventh-, and eighth-order Runge-Kutta formulas with stepsize control".
# National Aeronautics and Space Administration.
const fb_coefficients_78 = (8,
                            [     0      0      0       0        0         0       0       0     0      0    0 0
                                  2/27   0      0       0        0         0       0       0     0      0    0 0
                                  1/36   1/12   0       0        0         0       0       0     0      0    0 0
                                  1/24   0      1/8     0        0         0       0       0     0      0    0 0
                                  5/12   0    -25/16   25/16     0         0       0       0     0      0    0 0
                                  1/20   0      0       1/4      1/5       0       0       0     0      0    0 0
                                -25/108  0      0     125/108  -65/27    125/54    0       0     0      0    0 0
                                 31/300  0      0       0       61/225    -2/9    13/900   0     0      0    0 0
                                  2      0      0     -53/6    704/45   -107/9    67/90    3     0      0    0 0
                                -91/108  0      0      23/108 -976/135   311/54  -19/60   17/6  -1/12   0    0 0
                               2383/4100 0      0    -341/164 4496/1025 -301/82 2133/4100 45/82 45/164 18/41 0 0
                                  3/205  0      0       0        0        -6/41   -3/205  -3/41  3/41   6/41 0 0
                              -1777/4100 0      0    -341/164 4496/1025 -289/82 2193/4100 51/82 33/164 12/41 0 1],
                            # 7th order b-coefficients
                            [41/840 0 0 0 0 34/105 9/35 9/35 9/280 9/280 41/840 0 0],
                            # 8th order b-coefficients
                            [0 0 0 0 0 34/105 9/35 9/35 9/280 9/280 0 41/840 41/840],
                            )
ode78_fb(F, x0, tspan; kwargs...) = oderkf(F, x0, tspan, fb_coefficients_78...; kwargs...)

# Use Fehlberg version of ode78 by default
const ode78 = ode78_fb

# Use Dormand-Prince version of ode45 by default
const ode45 = ode45_dp

# Use Bogacki–Shampine version of ode23 by default
const ode23 = ode23_bs


# more higher-order embedded methods can be found in:
# P.J. Prince and J.R.Dormand: High order embedded Runge-Kutta formulae, Journal of Computational and Applied Mathematics 7(1), 1981.


#ODE4  Solve non-stiff differential equations, fourth order
#   fixed-step Runge-Kutta method.
#
#   [T,X] = ODE4(ODEFUN, X0, TSPAN) with TSPAN = [T0:H:TFINAL]
#   integrates the system of differential equations x' = f(t,x) from time
#   T0 to TFINAL in steps of H with initial conditions X0. Function
#   ODEFUN(T,X) must return a column vector corresponding to f(t,x). Each
#   row in the solution array X corresponds to a time returned in the
#   column vector T.
function ode4(F, x0, tspan)
    h = diff(tspan)
    x = Array(typeof(x0), length(tspan))
    x[1] = x0

    midxdot = Array(typeof(x0), 4)
    for i = 1:length(tspan)-1
        # Compute midstep derivatives
        midxdot[1] = F(tspan[i],         x[i])
        midxdot[2] = 2*F(tspan[i]+h[i]./2, x[i] + midxdot[1].*h[i]./2)
        midxdot[3] = 2*F(tspan[i]+h[i]./2, x[i] + midxdot[2].*h[i]./2)
        midxdot[4] = F(tspan[i]+h[i],    x[i] + midxdot[3].*h[i])

        # Integrate
        x[i+1] = x[i] + 1/6 .*h[i].*sum(midxdot)
    end
    return vcat(tspan), x
end

###############################################################################
## STIFF SOLVERS
###############################################################################

# Crude forward finite differences estimator of Jacobian as fallback

# FIXME: This doesn't really work if x is anything but a Vector or a scalar
function fdjacobian(F, x::Number, t)
    ftx = F(t, x)

    # The 100 below is heuristic
    dx = (x .+ (x==0))./100
    dFdx = (F(t,x+dx)-ftx)./dx

    return dFdx
end

function fdjacobian(F, x, t)
    ftx = F(t, x)
    lx = max(length(x),1)
    dFdx = zeros(eltype(x), lx, lx)
    for j = 1:lx
        # The 100 below is heuristic
        dx = zeros(eltype(x), lx)
        dx[j] = (x[j] .+ (x[j]==0))./100
        dFdx[:,j] = (F(t,x+dx)-ftx)./dx[j]
    end
    return dFdx
end

# ODE23S  Solve stiff systems based on a modified Rosenbrock triple
# (also used by MATLAB's ODE23s); see Sec. 4.1 in
#
# [SR97] L.F. Shampine and M.W. Reichelt: "The MATLAB ODE Suite," SIAM Journal on Scientific Computing, Vol. 18, 1997, pp. 1–22
#
# supports keywords: points = :all | :specified (using dense output)
#                    jacobian = G(t,y)::Function | nothing (FD)
function ode23s(F, y0, tspan; reltol = 1.0e-5, abstol = 1.0e-8,
                                                jacobian=nothing,
                                                points=:all,
                                                norm=Base.norm,
                                                minstep=abs(tspan[end] - tspan[1])/1e18,
                                                maxstep=abs(tspan[end] - tspan[1])/2.5,
                                                initstep=0.)


    # select method for computing the Jacobian
    if typeof(jacobian) == Function
        jac = jacobian
    else
        # fallback finite-difference
        jac = (t, y)->fdjacobian(F, y, t)
    end

    # constants
    const d = 1/(2 + sqrt(2))
    const e32 = 6 + sqrt(2)


    # initialization
    t = tspan[1]
    tfinal = tspan[end]
    tdir = sign(tfinal - t)

    h = initstep
    if h == 0.
      # initial guess at a step size
      h, F0 = hinit(F, y0, t, tdir, 3, reltol, abstol)
    else
      F0 = F(t,y0)
    end
    h = tdir*min(h, maxstep)

    y = y0
    tout = Array(typeof(t), 1)
    tout[1] = t         # first output time
    yout = Array(typeof(y0), 1)
    yout[1] = copy(y)         # first output solution


    J = jac(t,y)    # get Jacobian of F wrt y

    while abs(t) < abs(tfinal) && minstep < abs(h)
        if abs(t-tfinal) < abs(h)
            h = tfinal - t
        end

        if size(J,1) == 1
            W = one(J) - h*d*J
        else
            # note: if there is a mass matrix M on the lhs of the ODE, i.e.,
            #   M * dy/dt = F(t,y)
            # we can simply replace eye(J) by M in the following expression
            # (see Sec. 5 in [SR97])

            W = lufact( eye(J) - h*d*J )
        end

        # approximate time-derivative of F
        T = h*d*(F(t + h/100, y) - F0)/(h/100)

        # modified Rosenbrock formula
        k1 = W\(F0 + T)
        F1 = F(t + 0.5*h, y + 0.5*h*k1)
        k2 = W\(F1 - k1) + k1
        ynew = y + h*k2
        F2 = F(t + h, ynew)
        k3 = W\(F2 - e32*(k2 - F1) - 2*(k1 - F0) + T )

        err = (h/6)*norm(k1 - 2*k2 + k3) # error estimate
        delta = max(reltol*max(norm(y),norm(ynew)), abstol) # allowable error

        # check if new solution is acceptable
        if  err <= delta

            if points==:specified || points==:all
                # only points in tspan are requested
                # -> find relevant points in (t,t+h]
                for toi in tspan[(tspan.>t) & (tspan.<=t+h)]
                    # rescale to (0,1]
                    s = (toi-t)/h

                    # use interpolation formula to get solutions at t=toi
                    push!(tout, toi)
                    push!(yout, y + h*( k1*s*(1-s)/(1-2*d) + k2*s*(s-2*d)/(1-2*d)))
                end
            end
            if (points==:all) && (tout[end]!=t+h)
                # add the intermediate points
                push!(tout, t + h)
                push!(yout, ynew)
            end

            # update solution
            t = t + h
            y = ynew

            F0 = F2         # use FSAL property
            J = jac(t,y)    # get Jacobian of F wrt y
                            # for new solution
        end

        # update of the step size
        h = tdir*min( maxstep, abs(h)*0.8*(delta/err)^(1/3) )
    end

    return tout, yout
end


#ODEROSENBROCK Solve stiff differential equations, Rosenbrock method
#   with provided coefficients.
function oderosenbrock(F, x0, tspan, gamma, a, b, c; jacobian=nothing)

    if typeof(jacobian) == Function
        G = jacobian
    else
        G = (t, x)->fdjacobian(F, x, t)
    end

    h = diff(tspan)
    x = Array(typeof(x0), length(tspan))
    x[1] = x0

    solstep = 1
    while abs(tspan[solstep]) < abs(maximum(tspan))
        ts = tspan[solstep]
        hs = h[solstep]
        xs = x[solstep]
        dFdx = G(ts, xs)
        # FIXME
        if size(dFdx,1) == 1
            jac = 1/gamma/hs - dFdx[1]
        else
            jac = eye(dFdx)/gamma/hs - dFdx
        end

        g = Array(typeof(x0), size(a,1))
        g[1] = (jac \ F(ts + b[1]*hs, xs))
        x[solstep+1] = x[solstep] + b[1]*g[1]

        for i = 2:size(a,1)
            dx = zero(x0)
            dF = zero(x0/hs)
            for j = 1:i-1
                dx += a[i,j]*g[j]
                dF += c[i,j]*g[j]
            end
            g[i] = (jac \ (F(ts + b[i]*hs, xs + dx) + dF/hs))
            x[solstep+1] += b[i]*g[i]
        end
        solstep += 1
    end
    return vcat(tspan), x
end


# Kaps-Rentrop coefficients
const kr4_coefficients = (0.231,
                          [0              0             0 0
                           2              0             0 0
                           4.452470820736 4.16352878860 0 0
                           4.452470820736 4.16352878860 0 0],
                          [3.95750374663  4.62489238836 0.617477263873 1.28261294568],
                          [ 0               0                0        0
                           -5.07167533877   0                0        0
                            6.02015272865   0.1597500684673  0        0
                           -1.856343618677 -8.50538085819   -2.08407513602 0],)

ode4s_kr(F, x0, tspan; jacobian=nothing) = oderosenbrock(F, x0, tspan, kr4_coefficients...; jacobian=jacobian)

# Shampine coefficients
const s4_coefficients = (0.5,
                         [ 0    0    0 0
                           2    0    0 0
                          48/25 6/25 0 0
                          48/25 6/25 0 0],
                         [19/9 1/2 25/108 125/108],
                         [   0       0      0   0
                            -8       0      0   0
                           372/25   12/5    0   0
                          -112/125 -54/125 -2/5 0],)

ode4s_s(F, x0, tspan; jacobian=nothing) = oderosenbrock(F, x0, tspan, s4_coefficients...; jacobian=jacobian)

# Use Shampine coefficients by default (matching Numerical Recipes)
const ode4s = ode4s_s

const ms_coefficients4 = [ 1      0      0     0
                          -1/2    3/2    0     0
                          5/12  -4/3  23/12 0
                          -9/24   37/24 -59/24 55/24]

# ODE_MS Fixed-step, fixed-order multi-step numerical method
#   with Adams-Bashforth-Moulton coefficients
function ode_ms(F, x0, tspan, order::Integer)
    h = diff(tspan)
    x = Array(typeof(x0), length(tspan))
    x[1] = x0

    if 1 <= order <= 4
        b = ms_coefficients4
    else
        b = zeros(order, order)
        b[1:4, 1:4] = ms_coefficients4
        for s = 5:order
            for j = 0:(s - 1)
                # Assign in correct order for multiplication below
                #  (a factor depending on j and s) .* (an integral of a polynomial with -(0:s), except -j, as roots)
                p_int = polyint(poly(diagm(-[0:j - 1; j + 1:s - 1])))
                b[s, s - j] = ((-1)^j / factorial(j)
                               / factorial(s - 1 - j) * polyval(p_int, 1))
            end
        end
    end

    # TODO: use a better data structure here (should be an order-element circ buffer)
    xdot = similar(x)
    for i = 1:length(tspan)-1
        # Need to run the first several steps at reduced order
        steporder = min(i, order)
        xdot[i] = F(tspan[i], x[i])

        x[i+1] = x[i]
        for j = 1:steporder
            x[i+1] += h[i]*b[steporder, j]*xdot[i-(steporder-1) + (j-1)]
        end
    end
    return vcat(tspan), x
end

# Use order 4 by default
ode4ms(F, x0, tspan) = ode_ms(F, x0, tspan, 4)

end # module ODE
