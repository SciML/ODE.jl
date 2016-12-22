isdefined(Base, :__precompile__) && __precompile__()
# Ordinary Differential Equation Solvers

module ODE

using Polynomials
using Compat
using DiffEqBase

import DiffEqBase: solve

include("algorithm_types.jl")

## minimal function export list
# adaptive non-stiff:
export ode23, ode45, ode78
# non-adaptive non-stiff:
export ode4, ode4ms
# adaptive stiff:
export ode23s
# non-adaptive stiff:
export ode4s

# Common Interface
export ODEjlAlgorithm

## complete function export list: see runtests.jl

###############################################################################
## Coefficient Tableaus
###############################################################################

# Butcher Tableaus, or more generally coefficient tables
# see Hairer & Wanner 1992, p. 134, 166

abstract Tableau{Name, S, T<:Real}
# Name is the name of the tableau/method (a symbol)
# S is the number of stages (an int)
# T is the type of the coefficients
#
# TODO: have a type parameter which specifies adaptive vs non-adaptive
#
# For all types of tableaus it assumes fields:
# order::(Int...) # order of the method(s)
#
# For Runge-Kutta methods it assumes fields:
# a::Matrix{T}  # SxS matrix
# b::Matrix{T}  # 1 or 2 x S matrix (fixed step/ adaptive)
# c::Vector{T}  # S
#
# For a tableau:
#  c1  | a_11   ....   a_1s
#  .   | a_21 .          .
#  .   | a_31     .      .
#  .   | ....         .  .
#  c_s | a_s1  ....... a_ss
# -----+--------------------
#      | b_1     ...   b_s   this is the one used for stepping
#      | b'_1    ...   b'_s  this is the one used for error-checking

Base.eltype{N,S,T}(b::Tableau{N,S,T}) = T
order(b::Tableau) = b.order
# Subtypes need to define a convert method to convert to a different
# eltype with signature:
Base.convert{Tnew<:Real}(::Type{Tnew}, tab::Tableau) = error("Define convert method for concrete Tableau types")

###############################################################################
## HELPER FUNCTIONS
###############################################################################

# estimator for initial step based on book
# "Solving Ordinary Differential Equations I" by Hairer et al., p.169
function hinit{T}(F, x0, t0::T, tend, p, reltol, abstol)
    # Returns first step, direction of integration and F evaluated at t0
    tdir = sign(tend-t0)
    tdir==0 && error("Zero time span")
    tau = max(reltol*norm(x0, Inf), abstol)
    d0 = norm(x0, Inf)/tau
    f0 = F(t0, x0)
    d1 = norm(f0, Inf)/tau
    if d0 < 1e-5 || d1 < 1e-5
        h0 = 1e-6
    else
        h0 = 0.01*(d0/d1)
    end
    h0 = convert(T,h0)
    # perform Euler step
    x1 = x0 + tdir*h0*f0
    f1 = F(t0 + tdir*h0, x1)
    # estimate second derivative
    d2 = norm(f1 - f0, Inf)/(tau*h0)
    if max(d1, d2) <= 1e-15
        h1 = max(T(10)^(-6), T(10)^(-3)*h0)
    else
        pow = -(2 + log10(max(d1, d2)))/(p + 1)
        h1 = 10^pow
    end
    h1 = convert(T,h1)
    return tdir*min(100*h0, h1, tdir*(tend-t0)), tdir, f0
end

# isoutofdomain takes the state and returns true if state is outside
# of the allowed domain.  Used in adaptive step-control.
isoutofdomain(x) = isnan(x)

function make_consistent_types(fn, y0, tspan, btab::Tableau)
    # There are a few types involved in a call to a ODE solver which
    # somehow need to be consistent:
    #
    # Et = eltype(tspan)
    # Ey = eltype(y0)
    # Ef = eltype(Tf)
    #
    # There are also the types of the containers, but they are not
    # needed as `similar` is used to make containers.
    # Tt = typeof(tspan)
    # Ty = typeof(y0)              # note, this can be a scalar
    # Tf = typeof(F(tspan(1),y0))  # note, this can be a scalar
    #
    # Returns
    # - Et: eltype of time, needs to be a real "continuous" type, at
    #       the moment a AbstractFloat
    # - Eyf: suitable eltype of y and f(t,y)
    #   --> both of these are set to typeof(y0[1]/(tspan[end]-tspan[1]))
    # - Ty: container type of y0
    # - btab: tableau with entries converted to Et

    # Needed interface:
    # On components: /, -
    # On container: eltype, promote_type
    # On time container: eltype

    Ty = typeof(y0)
    Eyf = typeof(y0[1]/(tspan[end]-tspan[1]))

    Et = eltype(tspan)
    @assert Et<:Real
    if !(Et<:AbstractFloat)
        Et = promote_type(Et, Float64)
    end

    # if all are Floats, make them the same
    if Et<:AbstractFloat &&  Eyf<:AbstractFloat
        Et = promote_type(Et, Eyf)
        Eyf = Et
    end

    !isleaftype(Et) && warn("The eltype(tspan) is not a concrete type!  Change type of tspan for better performance.")
    !isleaftype(Eyf) && warn("The eltype(y0/tspan[1]) is not a concrete type!  Change type of y0 and/or tspan for better performance.")

    btab_ = convert(Et, btab)
    return Et, Eyf, Ty, btab_
end

###############################################################################
## NON-STIFF SOLVERS
###############################################################################

include("runge_kutta.jl")

# ODE_MS Fixed-step, fixed-order multi-step numerical method
#   with Adams-Bashforth-Moulton coefficients
function ode_ms(F, x0, tspan, order::Integer; kwargs...)
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
ode4ms(F, x0, tspan; kwargs...) = ode_ms(F, x0, tspan, 4; kwargs...)
ode5ms(F, x0, tspan; kwargs...) = ODE.ode_ms(F, x0, tspan, 5; kwargs...)

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

    h = initstep
    if h == 0.
        # initial guess at a step size
        h, tdir, F0 = hinit(F, y0, t, tfinal, 3, reltol, abstol)
    else
        tdir = sign(tfinal - t)
        F0 = F(t,y0)
    end
    h = tdir * min(abs(h), maxstep)

    y = y0
    tout = Array(typeof(t), 1)
    tout[1] = t         # first output time
    yout = Array(typeof(y0), 1)
    yout[1] = deepcopy(y)         # first output solution


    J = jac(t,y)    # get Jacobian of F wrt y

    while abs(t - tfinal) > 0 && minstep < abs(h)
        if abs(t-tfinal) < abs(h)
            h = tfinal - t
        end

        if size(J,1) == 1
            W = I - h*d*J
        else
            # note: if there is a mass matrix M on the lhs of the ODE, i.e.,
            #   M * dy/dt = F(t,y)
            # we can simply replace eye(J) by M in the following expression
            # (see Sec. 5 in [SR97])

            W = lufact( I - h*d*J )
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

        err = (abs(h)/6)*norm(k1 - 2*k2 + k3) # error estimate
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
function oderosenbrock(F, x0, tspan, gamma, a, b, c; jacobian=nothing, kwargs...)

    if typeof(jacobian) == Function
        G = jacobian
    else
        G = (t, x)->fdjacobian(F, x, t)
    end

    h = diff(tspan)
    x = Array(typeof(x0), length(tspan))
    x[1] = x0

    solstep = 1
    while solstep < length(tspan)
        ts = tspan[solstep]
        hs = h[solstep]
        xs = x[solstep]
        dFdx = G(ts, xs)

        jac = I/(gamma*hs) - dFdx

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

ode4s_kr(F, x0, tspan; jacobian=nothing, kwargs...) = oderosenbrock(F, x0, tspan, kr4_coefficients...; jacobian=jacobian, kwargs...)

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

ode4s_s(F, x0, tspan; jacobian=nothing, kwargs...) =
      oderosenbrock(F, x0, tspan, s4_coefficients...; jacobian=jacobian, kwargs...)

# Use Shampine coefficients by default (matching Numerical Recipes)
ode4s(F, x0, tspan; jacobian=nothing, kwargs...) = ode4s_s(F, x0, tspan; jacobian=nothing, kwargs...)

const ms_coefficients4 = [ 1      0      0     0
                          -1/2    3/2    0     0
                          5/12  -4/3  23/12 0
                          -9/24   37/24 -59/24 55/24]

####### Common Interface Bindings

include("common.jl")

end # module ODE
