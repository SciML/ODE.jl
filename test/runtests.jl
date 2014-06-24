using ODE
using Base.Test

tol = 1e-2

solvers = [
    ODE.ode23,
    ODE.ode23_bs,

    ODE.ode4,
    ODE.ode45_dp,
    ODE.ode45_fb,
    ODE.ode45_ck,

    ODE.ode4ms,
    ODE.ode4s_s,
    ODE.ode4s_kr
    ODE.ode78_fb]

for solver in solvers
    println("using $solver")
    # dy
    # -- = 6 ==> y = 6t
    # dt
    t,y=solver((t,y)->6, 0., [0:.1:1])
    @test maximum(abs(y-6t)) < tol

    # dy
    # -- = 2t ==> y = t.^2
    # dt
    t,y=solver((t,y)->2t, 0., [0:.001:1])
    @test maximum(abs(y-t.^2)) < tol
    
    # dy
    # -- = y ==> y = y0*e.^t
    # dt
    t,y=solver((t,y)->y, 1., [0:.001:1])
    @test maximum(abs(y-e.^t)) < tol

    # dv       dw
    # -- = -w, -- = v ==> v = v0*cos(t) - w0*sin(t), w = w0*cos(t) + v0*sin(t)
    # dt       dt
    #
    # y = [v, w]
    t,y=solver((t,y)->[-y[2], y[1]], [1., 2.], [0:.001:2*pi])
    ys = hcat(y...).'   # convert Vector{Vector{Float}} to Matrix{Float}
    @test maximum(abs(ys-[cos(t)-2*sin(t) 2*cos(t)+sin(t)])) < tol
end

println("All looks OK")
