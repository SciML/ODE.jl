export ode

type ODESolver{T<:Number}
	# The solver function with which to produce a solution
	solver::Function
	# The function f in y' = f(t,y)
	F::Function
	# The solution domain
	tspan::AbstractVector
	# Initial value for the solver
	x0::Union(T,AbstractVector{T})

	solve::Function
	
end

ODESolver(solver::Function, F::Function, tspan::AbstractVector, x0::Number) = ODESolver(solver, F, tspan, x0, () -> solver(F, tspan, x0))
ODESolver{T}(solver::Function, F::Function, tspan::AbstractVector, x0::AbstractVector{T}) = ODESolver(solver, F, tspan, x0, () -> solver(F, tspan, x0))

function ode{T<:Number}(solver::Symbol, F::Function, tspan::AbstractVector, x0::Union(T,AbstractVector{T}))
	if !contains((:ode23, :ode4, :ode45, :ode4s, :ode4ms), solver)
		error("$solver is not one of the valid ODE solvers")
	end

	ode = ODESolver(eval(solver), F, tspan, x0)
end
