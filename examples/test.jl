include("../src/ODE.jl")

module Test

using ODE

T = Float64
Y = Vector{T}
t0 = zero(T)
y0 = T[one(T)]

steppers = [# ODE.RKStepperAdaptive{:rk45},
            # ODE.RKStepperFixed{:feuler},
            ODE.DenseStepper]

for st in steppers
    ode  = ODE.ExplicitODE(t0,y0,(t,y,dy)->dy[1]=y[1])
    opts = Dict(:initstep=>0.1,
                :tspan=>[0.,0.5,1.],
                :points=>:specified,
                :reltol=>1e-5,
                :abstol=>1e-5)

    sol = ODE.solve(ode,st;opts...)

    println("Raw iterator")
    for (t,y) in sol
        println((t,y,norm(y-[exp(t)])))
    end

    println(collect(sol))
end

end
