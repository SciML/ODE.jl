# Ordinary Differential Equation Solvers

"""
Coding conventions:

- use t::T, y::Y, dy::Y
- p::Problem, ::P
- ivp::IVP, ::O
- integrator::AbstractIntegrator, ::I
- opts::AbstactOptions

Variables and Type variables:
- T -> t::T
- Y -> y::Y  TODO: or Vector{Y}?


"""
module ODE

using Polynomials
using Compat
import Compat.String
using ForwardDiff

import Base: start, next, done, collect, show, convert

# Core infrastructure
#
# When wrapping a new solver it will need to use and conform to
# methods and types within these files.
#
# Note, if we go the MathProgBase.jl route, then these files would go
# into ODEBase.jl.
include("types.jl")
include("tableaus.jl")
include("options.jl")
include("helpers.jl")

# Dense output wrapper
include("dense.jl")

# Particular solvers
include("ode23s.jl")
include("runge-kutta.jl")
include("adams-bashford-moulton.jl")
include("rosenbrock.jl")

# User interface to solvers
include("interfaces.jl")

end # module ODE
