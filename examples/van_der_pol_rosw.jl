# Van der Pol using DASSL.jl
using ODE
using DASSL
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
    # d vdp_impl /dy + α d vdp_impl/dydot
    #
    # Jacobian
    return Jvdp(0,y) - α * eye(2)
end

# Results 24dfdc9cacfd6
# elapsed time: 0.025767095 seconds (25194952 bytes allocated)
# elapsed time: 0.00997568 seconds   (6674864 bytes allocated)
# Adaptive step: abs error of ode23s vs ref:
# [0.012513592025336528,0.013225223452835055]
# Adaptive step: abs error of ode_rosw vs ref:
# [0.012416936355840624,0.013124961519338618]


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
    # d vdp_impl /dy + α d vdp/dydot
    #
    # Jacobian

    res[1,1] = 0 - α
    res[1,2] = 1
    res[2,1] = -mu^2*(y[2]*2*y[1] + 1)
    res[2,2] = mu^2*(1-y[1]^2) - α
    return nothing
end

# elapsed time: 0.02611999 seconds (25194952 bytes allocated)
# elapsed time: 0.006630583 seconds (4025352 bytes allocated)
# elapsed time: 0.005866466 seconds (2670408 bytes allocated) (newer version)
# Adaptive step: abs error of ode23s vs ref:
# [0.012513592025336528,0.013225223452835055]
# Adaptive step: abs error of ode_rosw vs ref:
# [0.012416936355840624,0.013124961519338618]


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

t,yout3 = ode23s(vdp, y0, [0,0.1]; jacobian=Jvdp)
gc()
@time t,yout3 = ode23s(vdp, y0, tspan; jacobian=Jvdp)

t,yout4 = ode_rosw(vdp_impl!, Jvdp_impl!, y0, [0, 0.1])
gc()
@time t,yout4 = ode_rosw(vdp_impl!, Jvdp_impl!, y0, tspan)

println("Adaptive step: abs error of ode23s vs ref:")
println(yout3[end]-refsol)
println("Adaptive step: abs error of ode_rosw vs ref:")
println(yout4[end]-refsol)


####################
# DAE: reduced Van der Pol



# implicit inplace equation
function dae_vdp_impl!(res, y, ydot)
    #  dae_vdp_impl(t, y)
    #
    # Implicit, reduced van der Pol equation.  mu->\infinity

    res[1] = y[2] - ydot[1]
    res[2] = y[1] - (y[2]^3/3 - y[2])
    return nothing
end
function Jdae_vdp_impl!(res, y, ydot, α)
    # d dae_vdp_impl /dy + α d dae_vdp_impl/dydot
    #
    # Jacobian

    res[1,1] = 0 - α
    res[1,2] = 1
    res[2,1] = 1
    res[2,2] = -(y[2]^2 - 1)
    return nothing
end

t,yout5 = ode_rosw(dae_vdp_impl!, Jdae_vdp_impl!, y0, [0, 0.1])
gc()
@time t,yout5 = ode_rosw(dae_vdp_impl!, Jdae_vdp_impl!, y0, tspan)



# implicit inplace equation
function dae_vdp(t, y, ydot)
    #  dae_vdp_impl(t, y)
    #
    # Implicit, reduced van der Pol equation.  mu->\infinity

    [y[2] - ydot[1], y[1] - (y[2]^3/3 - y[2])]
end
function Jdae_vdp(t, y, ydot, α)
    # d dae_vdp_impl /dy + α d dae_vdp_impl/dydot
    #
    # Jacobian
    res = zeros(2,2)
    res[1,1] = 0 - α
    res[1,2] = 1
    res[2,1] = 1
    res[2,2] = -(y[2]^2 - 1)
    return res
end
tspan = [0, 4e6]
y0 = [1., 0, 0]

t,ydassl = dasslSolve(dae_vdp, y0, tspan)



