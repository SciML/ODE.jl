using ODE, DiffEqBase, DiffEqProblemLibrary, Test
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear
using LinearAlgebra

dt=1/2^(4)

algs = [ode23(),ode45(),ode78(),ode4(),ode4ms(),ode4s()] # no ode23s

# Check for errors

prob = prob_ode_linear

for alg in algs
    sol =solve(prob,alg;dt=dt,abstol=1e-6,reltol=1e-3)
    @test typeof(sol[2]) <: Number
end

prob = prob_ode_2Dlinear

for alg in algs
    sol =solve(prob,alg;dt=dt,dtmin=eps(),abstol=1e-6,reltol=1e-3)
    @test size(sol[2]) == (4,2)
end

# Check that warnings are displayed
for alg in algs
    sol = solve(prob,alg;dt=dt,abstol=1e-6,reltol=1e-3, dense=true, verbose=true,
                save_idxs = true, progress = true, beta1 = 1.23, gamma = 0.5)
end
