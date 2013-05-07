
module DETEST

export IVPWithSolution, single_equations

# The test problems and solutions are taken from from http://dx.doi.org/10.1137/0709052

abstract AbstractIVP{T<:Union(Number, AbstractVector{Number})}

type IVPWithSolution{T<:Union(Float64, AbstractVector{Float64})} <: AbstractIVP
    f::Function                     # ODE definition    y' = f(t,y)
    y0::T                           # Initial condition y(0) = 0 
    tspan::AbstractVector{Float64}   # Boundaries for time domain
    sol::Function                   # Analytic solution y(t)
end

# Single equations
single_equations = [
    # A1: The negative potential
    IVPWithSolution( (t,y)->-y , 1. , [0., 10.] , t->exp(-t) )
    # A2: A special case of the Riccati equation
#    IVPWithSolution((t,y)->-y^3 / 2, 1, [0, 5], t-> 1/sqrt(t+1))
]



end