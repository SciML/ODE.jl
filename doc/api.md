#General API for all solvers
Please open pull requests or issues to propose changes or clarifications.

##Basic interface
The general interface for all ODE solvers is:

```julia
t_out, y_out = odeXX(F, y0, tspan; kwargs...)
```

Each (adaptive) solver accepts 3 arguments

- `F`: the RHS of the ODE `dy/dt = F(t,y)`, which is a function of `t` and `y(t)` and returns `dy/dt::typeof(y./t)`
- `y0`: initial value for `y`. The type of `y0` determines the element type of the `y_out` vector (`y_out::Vector{typeof(y0)}`)
- `tspan`: Any iterable of sorted `t` values at which the solution (`y`) is requested. Most solvers will only consider `tspan[1]` and `tspan[end]`, and intermediary points will be interpolated. If `tspan[1] > tspan[end]` the integration is performed backwards.

The solver returns two arrays

- `tout`: Vector of points at which solutions were obtained (also see keyword `points`)
- `yout`: solutions at times `tout`; if `points=:all | :specified` is used, `yout::Vector{typeof(y0)}`. Note that if `y0` is a vector, you can get a matlab-like matrix with `hcat(yout...)`.

Each solver might implement its own keywords, but the following keywords have a standardized interpretation across all solvers. Solvers should raise an error if a unrecognized keyword argument is used.

- `norm`: user-supplied norm for determining the error `E` (default `Base.norm`)
- `abstol` and/or `reltol`: an integration step is accepted if `E <= abstol || E <= reltol*abs(y)` (ideally we want both criteria for all solvers, **done** in #13)
- `points=:all | :specified`: controls the type of output according to
 * `points==:all` (default) output is given for each value in `tspan` as well as for each intermediate point the solver used. 
 * `points==:specified` output is given only for each value in `tspan`.
- `maxstep`, `minstep` and `initstep`: determine the maximal, minimal and initial integration step.
- `retries = 0` Sometimes an integration step takes you out of the region where `F(t,y)` has a valid solution and F might throw `DomainError` or other exceptions. `retries` sets a limit to the number of times the solver might try with a smaller step. 

##Iterator interface
Under construction #27
