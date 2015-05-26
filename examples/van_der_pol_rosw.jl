# Van der Pol using DASSL.jl
using ODE
abstol = 1e-6
reltol = 1e-6

const mu = 100.

# Explicit equation
function vdp(t, y)
    #  vdp(t, y)
    #
    # ydot of van der Pol equation.  (consider using rescaled eq. so
    # period stays the same for different mu.  Cf. Wanner & Hairer 1991,
    # p.5)
    return [y[2], mu^2*((1-y[1]^2) * y[2] - y[1])]
end
function Jvdp(t, y)
    #  Jvdp(t, y, mu)
    #
    # Jacobian of vdp
    return [0 1;
            -mu^2*(y[2]*2*y[1] + 1)  mu^2*(1-y[1]^2)]
end

# implicit equation
function vdp_impl(y, ydot)
    #  vdp_impl(t, y)
    #
    # Implicit van der Pol equation.  (consider using rescaled eq. so
    # period stays the same for different mu.  Cf. Wanner & Hairer 1991,
    # p.5)
    return vdp(0,y) - ydot
end
function Jvdp_impl(y, ydot, α)
    # d vdp_impl /dy + α d vdp/dy
    #
    # Jacobian
    return Jvdp(0,y) - α * eye(2)
end

# implicit inplace equation
function vdp_impl!(res, y, ydot)
    #  vdp_impl(t, y)
    #
    # Implicit van der Pol equation.  (consider using rescaled eq. so
    # period stays the same for different mu.  Cf. Wanner & Hairer 1991,
    # p.5)
    #
    # Note, that ydot can be used as res here.
    res[1] = y[2] - ydot[1]
    res[2] = mu^2*((1-y[1]^2) * y[2] - y[1]) - ydot[2]
    return nothing
end
function Jvdp_impl!(res, y, ydot, α)
    # d vdp_impl /dy + α d vdp/dy
    #
    # Jacobian

    res[1,1] = 0 - α
    res[1,2] = 1
    res[2,1] = -mu^2*(y[2]*2*y[1] + 1)
    res[2,2] = mu^2*(1-y[1]^2) - α
    return nothing
end


###
# reference as t=2
##
refsol = [0.1706167732170483e1, -0.8928097010247975e0]
y0 = [2.,0.];
tstart = 0.
tend = 2.



###
# Fixed step
###
# nt = 50_000;
# tspan = linspace(tstart, tend, nt)

# t,yout1 = ode4s(vdp, y0, tspan; jacobian=Jvdp)
# @time t,yout1 = ode4s(vdp, y0, tspan; jacobian=Jvdp)

# t,yout2 = ode_rosw_fixed(vdp_impl, Jvdp_impl, y0, tspan)
# @time t,yout2 = ode_rosw_fixed(vdp_impl, Jvdp_impl, y0, tspan)

# println("Fixed step: abs error of ode4s vs ref:")
# println(yout1[end]-refsol)
# println("Fixed step: abs error of ode_rosw vs ref:")
# println(yout2[end]-refsol)

# #### adaptive
tspan = linspace(tstart, tend, 2)

t,yout3 = ode23s(vdp, y0, tspan; jacobian=Jvdp)
gc()
@time t,yout3 = ode23s(vdp, y0, tspan; jacobian=Jvdp)

t,yout4 = ode_rosw(vdp_impl, Jvdp_impl, y0, tspan)
gc()
@time t,yout4 = ode_rosw(vdp_impl, Jvdp_impl, y0, tspan)

println("Adaptive step: abs error of ode23s vs ref:")
println(yout3[end]-refsol)
println("Adaptive step: abs error of ode_rosw vs ref:")
println(yout4[end]-refsol)

# rosw
# @time yout3, ts, steps, dts, xerrs = rosw_runner_adaptive(
#                      vdp_impl, Jvdp_impl, [0., 2.], y0;
#                      reltol=reltol, abstol=abstol, dt0=1e-5, mindt=1e-8)

# check

# using Winston
# plot(ts, yout3[1,:])

# oplot(tspan, yout2[1,:],"r")
