include("../src/ODE.jl")

module Test

using ODE

T = Float64
Y = Vector{T}
t0 = zero(T)
y0 = T[one(T)]

steppers = [ODE.RKStepperAdaptive{:rk45},
            ODE.RKStepperFixed{:feuler}]

for st in steppers
    ode  = ODE.ExplicitODE(t0,y0,(t,y,dy)->dy[1]=y[1])
    opts = Dict(:initstep=>0.1,
                :tstop=>1.,
                # :tspan=>[0.,1.],
                :points=>:specified,
                :reltol=>1e-5,
                :abstol=>1e-5)

    sol = ODE.solve(ode,st;opts...)
    den = ODE.dense(sol;opts...)
    println(typeof(sol.stepper.options))
    println(sol.stepper.options)

    println("Raw iterator")
    for (t,y) in sol
        println((t,y))
    end

    println("Dense output")
    for (t,y) in den
        println((t,y))
    end

    println(collect(sol))
    println(collect(den))
end

end
