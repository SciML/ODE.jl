using ODE
using TestProblems
using FactCheck

abstol = 1e-6
reltol = 1e-6

#require("nearequal")

#println("Running tests for package ODE")

problem_set = [
    TestProblems.ODETest(TestProblems.A1, ode45, 1e-6, 1e-6)
]

@facts begin

    for problem in problem_set
        tout, yout = problem.solver(problem.ivp.f, problem.ivp.tspan, problem.ivp.y0)
        soly = reshape(map(problem.ivp.sol,tout),size(yout))

        @fact yout => roughly(soly, reltol)
    end
end
