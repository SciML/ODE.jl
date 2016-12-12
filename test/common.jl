using ODE, DiffEqBase, DiffEqProblemLibrary

dt=1/2^(4)

prob = prob_ode_linear

sol =solve(prob,ODE.ode23Alg();dt=dt)
# plot(sol,plot_analytic=true)

sol =solve(prob,ode23sAlg();dt=dt)
sol =solve(prob,ode45Alg();dt=dt)
sol =solve(prob,ode78Alg();dt=dt)

prob = prob_ode_2Dlinear

sol =solve(prob,ode23Alg();dt=dt,dtmin=eps())
# plot(sol,plot_analytic=true)

sol =solve(prob,ode23sAlg();dt=dt)
sol =solve(prob,ode45Alg();dt=dt)
sol =solve(prob,ode78Alg();dt=dt)
