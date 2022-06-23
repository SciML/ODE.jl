# Here are tests which test what interface the solvers require.

################################################################################
# This is to test a scalar-like state variable
# (due to @acroy: https://gist.github.com/acroy/28be4f2384d01f38e577)

import Base: +, -, *, /
using LinearAlgebra

const delta0 = 0.0
const V0 = 1.0
const g0 = 0.0

# define custom type ...
struct CompSol <: AbstractMatrix{ComplexF64}
    rho::Matrix{ComplexF64}
    x::Float64
    p::Float64

    CompSol(r, x, p) = new(copy(r), x, p)
end

Base.getindex(A::CompSol, i...) = A.rho[i...]
Base.setindex!(A, x, i...) = setindex!(A.rho, x, i...)

# ... which has to support the following operations
# to work with odeX
function LinearAlgebra.norm(y::CompSol, p::Float64)
    maximum([LinearAlgebra.norm(y.rho, p) abs(y.x) abs(y.p)])
end
LinearAlgebra.norm(y::CompSol) = LinearAlgebra.norm(y::CompSol, 2.0)

+(y1::CompSol, y2::CompSol) = CompSol(y1.rho + y2.rho, y1.x + y2.x, y1.p + y2.p)
-(y1::CompSol, y2::CompSol) = CompSol(y1.rho - y2.rho, y1.x - y2.x, y1.p - y2.p)
*(y1::CompSol, s::Real) = CompSol(y1.rho * s, y1.x * s, y1.p * s)
*(s::Real, y1::CompSol) = y1 * s
/(y1::CompSol, s::Real) = CompSol(y1.rho / s, y1.x / s, y1.p / s)

### new for PR #68
Base.abs(y::CompSol) = LinearAlgebra.norm(y, 2.0) # TODO not needed anymore once https://github.com/JuliaLang/julia/pull/11043 is in current stable julia
Base.zero(::Type{CompSol}) = CompSol(complex(zeros(2, 2)), 0.0, 0.0)
ODE.isoutofdomain(y::CompSol) = any(isnan, vcat(y.rho[:], y.x, y.p))

################################################################################

# define RHSs of differential equations
# delta, V and g are parameters
function rhs(t, y, delta, V, g)
    H = [[-delta / 2 V]; [V delta / 2]]

    rho_dot = -im * H * y.rho + im * y.rho * H
    x_dot = y.p
    p_dot = -y.x

    return CompSol(rho_dot, x_dot, p_dot)
end

# inital conditons
rho0 = zeros(2, 2);
rho0[1, 1] = 1.0;
y0 = CompSol(complex(rho0), 2.0, 1.0)

# solve ODEs
endt = 2 * pi;

t, y1 = ODE.ode45((t, y) -> rhs(t, y, delta0, V0, g0), y0, [0.0, endt]) # used as reference
print("Testing interface for scalar-like state... ")
for solver in solvers
    # these only work with some Array-like interface defined:
    if solver in [ODE.ode23s, ODE.ode4s_s, ODE.ode4s_kr]
        continue
    end
    t, y2 = solver((t, y) -> rhs(t, y, delta0, V0, g0), y0,
                   range(0.0, stop = endt, length = 500))
    @test LinearAlgebra.norm(y1[end] - y2[end]) < 0.15
end
println("ok.")

################################################################################
# TODO: test a vector-like state variable, i.e. one which can be indexed.
