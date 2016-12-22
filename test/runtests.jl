using ODE
using Base.Test

tols = [5e-2, 1e-2, 1e-2]

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
           # adaptive
           ODE.ode23s]

# Because of https://github.com/JuliaLang/julia/issues/16667
# which was fixed for BigFloat in https://github.com/JuliaLang/julia/pull/16999
# but not back-ported to Julia 0.4
if VERSION<=v"0.5-"
    dtypes = [Float32, Float64]
else
    dtypes = [Float32, Float64, BigFloat]
end
for solver in solvers
    println("using $solver")
    for (tol,T) in zip(tols, dtypes)
        if solver==ODE.ode45_fe && T==Float32
            # For some reason ode45_fe hangs with Float32!
            continue
        end

        # dy
        # -- = 6 ==> y = 6t
        # dt
        t,y=solver((t,y)->T(6), T(0), T[0:.1:1;])
        @test maximum(abs(y-6t)) < tol
        @test eltype(t)==T

        # dy
        # -- = 2t ==> y = t.^2
        # dt
        t,y=solver((t,y)->T(2)*t, T(0), T[0:.001:1;])
        @test maximum(abs(y-t.^2)) < tol
        @test eltype(t)==T

        # dy
        # -- = y ==> y = y0*e.^t
        # dt
        t,y=solver((t,y)->y, T(1), T[0:.001:1;])
        @test maximum(abs(y-e.^t)) < tol
        @test eltype(t)==T

        t,y=solver((t,y)->y, T(1), T[1:-.001:0;])
        @test maximum(abs(y-e.^(t-1))) < tol
        @test eltype(t)==T

        # dv       dw
        # -- = -w, -- = v ==> v = v0*cos(t) - w0*sin(t), w = w0*cos(t) + v0*sin(t)
        # dt       dt
        #
        # y = [v, w]
        t,y=solver((t,y)->[-y[2]; y[1]], T[1, 2], T[0:.001:2*pi;])
        ys = hcat(y...).'   # convert Vector{Vector{T}} to Matrix{T}
        @test maximum(abs(ys-[cos(t)-2*sin(t) 2*cos(t)+sin(t)])) < tol
        @test eltype(t)==T
    end
end

# Test negative starting times ODE.ode23s
@test length(ODE.ode23s((t,y)->[-y[2]; y[1]], [1., 2.], [-5., 0])[1]) > 1


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
include("common.jl")
println("All looks OK")
