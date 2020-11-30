Various basic Ordinary Differential Equation solvers implemented in Julia.

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://github.com/SciML/ODE.jl/workflows/CI/badge.svg)](https://github.com/SciML/ODE.jl/actions?query=workflow%3ACI)
[![Coverage Status](https://img.shields.io/coveralls/JuliaDiffEq/ODE.jl.svg)](https://coveralls.io/r/JuliaDiffEq/ODE.jl)

Pull requests are always highly welcome to fix bugs, add solvers, or anything else!

# Current status of the project

This project is deprecated in favor of [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl) and its ODE solvers [OrdinaryDiffEq.jl](https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl). This library is in "maitanance mode", meaning that it is being upgraded with each Julia version, but not seeing active feature development. ODE.jl contains the basic functionality that was moved here when the package was originally moved from Base. Although quite poorly tested, at least some of the functionality is quite reliable. Use at your own risk.

## Usage On the Common Interface

The ODE.jl methods can be used on the common interface. Simply use the solver's name as the algorithm. For example, [the ODE tutorial](http://docs.juliadiffeq.org/dev/tutorials/ode_example.html) can be solved using ODE.jl's `ode45` by using the following commands:

```julia
using ODE
f(u,p,t) = 1.01*u
u0=1/2
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob,ode45(),reltol=1e-8,abstol=1e-8)
using Plots
plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
     xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!") # legend=false
plot!(sol.t, t->0.5*exp(1.01t),lw=3,ls=:dash,label="True Solution!")
```

Note that ODE.jl does not natively support inplace updates. Inplace functions `f(t,u,du)` are converted to out-of-place functions `du=f(t,u)` and thus it will not be any more efficient.

## Basic API

All of the ODE.jl solvers the following basic API:

    tout, yout = odeXX(F, y0, tspan; keywords...)

to solve the explicitly defined ODE by dy/dt = F(t,y). A few other solvers are also exported, see the source code for details.

The adaptive solvers accept the following keywords
- `norm`: user-supplied norm for determining the error `E` (default `Base.vecnorm`),
- `abstol` and/or `reltol`: an integration step is accepted if `E <= abstol || E <= reltol*abs(y)` (defaults `reltol = 1e-5`, `abstol = 1e-8`),
- `maxstep`, `minstep` and `initstep`: determine the maximal, minimal and initial integration step (defaults `minstep=|tspan[end] - tspan[1]|/1e9`, `maxstep=|tspan[end] - tspan[1]|/2.5` and automatic initial step estimation).
- `points=:all` (default): output is given for each value in `tspan` as well as for each intermediate point the solver used.
- `points=:specified`: output is given only for each value in `tspan`.

Additionally, `ode23s` solver supports
- `jacobian = G(t,y)`: user-supplied Jacobian G(t,y) = dF(t,y)/dy (default estimate by finite-difference method).

There are also fixed step Runge-Kutta and Rosenbrock solvers available.

## Available Solvers

Currently, `ODE` exports the following adaptive solvers:

* `ode23`: 2nd order adaptive solver with 3rd order error control, using the Bogacki–Shampine coefficients
* `ode45`: 4th order adaptive solver with 5th order error control, using the Dormand Prince coefficients. Fehlberg and Cash-Karp coefficients are also available.
* `ode78`: 7th order adaptive solver with 8th order error control, using the Fehlberg coefficients.
* `ode23s`: 2nd/3rd order adaptive solver for stiff problems, using a modified Rosenbrock triple.

For a full list, see the [DiffEqDocs ODE Solvers page](http://docs.juliadiffeq.org/dev/solvers/ode_solve.html#ODE.jl-1).

# Examples
The examples directory contain a few notebooks that show how to get started. You can also see them here:
* [Simple differential equation](http://nbviewer.jupyter.org/github/JuliaLang/ODE.jl/blob/master/examples/Simple_Differential_Equation.ipynb)
* [Lorenz Attractor](http://nbviewer.jupyter.org/github/JuliaLang/ODE.jl/blob/master/examples/Lorenz_Attractor.ipynb)
* [Terminal Velocity](http://nbviewer.jupyter.org/github/JuliaLang/ODE.jl/blob/master/examples/Terminal_Velocity.ipynb)
