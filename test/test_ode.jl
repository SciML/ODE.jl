using ODE
using Test

tol = 1e-2

solvers = [
    ODE.ode23,
    ODE.ode4,
    ODE.ode45_dp,
    ODE.ode45_fb,
    ODE.ode45_ck,
    
    ODE.ode4ms,
    ODE.ode4s_s,
    ODE.ode4s_kr]

for solver in solvers
    println("using $solver")
    # dy
    # -- = 6 ==> y = 6t
    # dt
    t,y=solver((t,y)->6, [0:.1:1], [0.])
    @test max(abs(y-6t)) < tol

    # dy
    # -- = 2t ==> y = t.^2
    # dt
    t,y=solver((t,y)->2t, [0:.001:1], [0.])
    @test max(abs(y-t.^2)) < tol
    
    # dy
    # -- = y ==> y = y0*e.^t
    # dt
    t,y=solver((t,y)->y, [0:.001:1], [1.])
    @test max(abs(y-e.^t)) < tol
end



println("All looks OK")
