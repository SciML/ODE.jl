include("ODE.jl")

module Test

using ODE

T = Float64
Y = Vector{T}
t0 = zero(T)
y0 = T[one(T)]

st   = ODE.RKStepperAdaptive{:rk45}
ode  = ODE.ExplicitODE(t0,y0,(t,y,dy)->dy[1]=y[1])
opts = Dict(:initstep=>0.1,
            :tstop=>1.,
            :tspan=>[0.,1.],
            :points=>:specified,
            :reltol=>1e-5,
            :abstol=>1e-5)

stepper = st{T}(ode)
sol = ODE.Solver(ode,stepper)
println(sol)


sol = ODE.solve(ode,st;opts...)
den = ODE.dense(sol;opts...)
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

println(collect(sol'))
println(collect(den'))

end
