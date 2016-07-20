# Ordinary Differential Equation Solvers

module ODE

using Polynomials
using Compat
using ForwardDiff

import Base.convert, Base.show
import Base: start, next, done, call, collect

# basic type definitions
include("types.jl")
include("options.jl")
include("helpers.jl")

# dense output wrapper
include("dense.jl")

# particular solvers
include("ode23s.jl")
include("runge-kutta.jl")
# include("multistep.jl")
include("adams-bashford-moulton.jl")
include("rosenbrock.jl")
# include("taylor.jl")

include("iterators.jl")
include("interfaces.jl")

end # module ODE
