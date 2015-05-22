using ODE
using Base.Test

tol = 1e-2

solvers = [
           ## Non-stiff           
           # fixed step
           ODE.ode1,
           ODE.ode2_midpoint,
           ODE.ode2_heun,
           ODE.ode4,
           ODE.ode4ms,
           ODE.ode5ms,
           # adaptive
#           ODE.ode21, # this fails on Travis with 0.4?! TODO revert once fixed.
           ODE.ode23,
           ODE.ode45_dp,
           ODE.ode45_fe,
           ODE.ode78,

           ## Stiff
           # fixed-step
           ODE.ode4s_s,
           ODE.ode4s_kr,
#           ODE.ode_rosw,
           # adaptive
           ODE.ode23s]

for solver in solvers
    println("using $solver")
    # dy
    # -- = 6 ==> y = 6t
    # dt
    t,y=solver((t,y)->6.0, 0., [0:.1:1;])
    @test maximum(abs(y-6t)) < tol

    # dy
    # -- = 2t ==> y = t.^2
    # dt
    t,y=solver((t,y)->2t, 0., [0:.001:1;])
    @test maximum(abs(y-t.^2)) < tol
    

    # dy
    # -- = y ==> y = y0*e.^t
    # dt
    t,y=solver((t,y)->y, 1., [0:.001:1;])
    @test maximum(abs(y-e.^t)) < tol

    t,y=solver((t,y)->y, 1., [1:-.001:0;])
    @test maximum(abs(y-e.^(t-1))) < tol
 
    # dv       dw
    # -- = -w, -- = v ==> v = v0*cos(t) - w0*sin(t), w = w0*cos(t) + v0*sin(t)
    # dt       dt
    #
    # y = [v, w]
    t,y=solver((t,y)->[-y[2]; y[1]], [1., 2.], [0:.001:2*pi;])
    ys = hcat(y...).'   # convert Vector{Vector{Float}} to Matrix{Float}
    @test maximum(abs(ys-[cos(t)-2*sin(t) 2*cos(t)+sin(t)])) < tol
end

# rober testcase from http://www.unige.ch/~hairer/testset/testset.html
let
    println("ROBER test case")
    function f(t, y)
        ydot = similar(y)
        ydot[1] = -0.04*y[1] + 1.0e4*y[2]*y[3]
        ydot[3] = 3.0e7*y[2]*y[2]
        ydot[2] = -ydot[1] - ydot[3]
        ydot
    end
    t = [0., 1e11]
    t,y = ode23s(f, [1.0, 0.0, 0.0], t; abstol=1e-8, reltol=1e-8,
                                        maxstep=1e11/10, minstep=1e11/1e18)

    refsol = [0.2083340149701255e-07,
              0.8333360770334713e-13,
              0.9999999791665050] # reference solution at tspan[2]
    @test norm(refsol-y[end], Inf) < 2e-10
end
include("interface-tests.jl")

println("All looks OK")
