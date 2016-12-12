using ODE, DiffEqBase, DiffEqProblemLibrary

dt=1/2^(4)

prob = prob_ode_linear

sol =solve(prob,ODE.ode23();dt=dt)
# plot(sol,plot_analytic=true)

sol =solve(prob,ode23s();dt=dt)
sol =solve(prob,ode45();dt=dt)
sol =solve(prob,ode78();dt=dt)

prob = prob_ode_2Dlinear

sol =solve(prob,ode23();dt=dt,dtmin=eps())
# plot(sol,plot_analytic=true)

sol =solve(prob,ode23s();dt=dt)
sol =solve(prob,ode45();dt=dt)
sol =solve(prob,ode78();dt=dt)
