using ODE

# Define IVP-instance which holds the mathematical problem definition:
t0 = 0.0
y0 = [1.0]
ode  = ODE.ExplicitODE(t0,y0,(t,y,dy)->dy[1]=y[1])

# options for the solver
opts = Dict(:initstep=>0.1,
            :tstop=>1.0,
            :reltol=>1e-5,
            :abstol=>1e-5)
# pick your solver
stepper = [ODE.RKStepperAdaptive{:rk45},
           ODE.ModifiedRosenbrockStepper][2]

# create a Solver instance
sol = ODE.solve(ode,stepper;opts...)

 # iterate over the solution
println("t, y, err")
for (t,y) in sol
    println((t,y[1],abs(y[1]-e.^t)))
end

# or collect it
println(collect(sol))

### Reverse time integration, rest as above
t0 = 1.0
y0 = [1.0]
ode  = ODE.ExplicitODE(t0,y0,(t,y,dy)->dy[1]=y[1])
opts = Dict(:initstep=>0.1,
            :tstop=>0.0,
            :reltol=>1e-5,
            :abstol=>1e-5)

sol = ODE.solve(ode,stepper;opts...)

println("t, y, err")
for (t,y) in sol # iterate over the solution
    println((t,y[1],abs(y[1]-e.^(t-1))))
end

println(collect(sol))
