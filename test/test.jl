require("../src/ODE.jl")
require("DETEST.jl")
#require("src/nearequal.jl")

using ODE
using DETEST
using FactCheck

abstol = 1e-6
reltol = 1e-6

#require("nearequal")

#println("Running tests for package ODE")

solvers = [ode4]
problem = single_equations[1]

@facts "ode45" begin

    tout,yout = ode45(problem.f, problem.tspan, problem.y0)
    soly = reshape(map(problem.sol,tout),size(yout))

    @fact yout => roughly(soly, reltol)
end
