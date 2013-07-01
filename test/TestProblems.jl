
module TestProblems

export IVPWithSolution, ODETest, A1, A2

type IVPWithSolution
    name::String
    f::Function                     # ODE definition    y' = f(t,y)
    y0::Number                      # Initial condition y(0) = 0 
    tspan::AbstractVector{Number}   # Boundaries for time domain
    sol::Function                   # Analytic solution y(t)
end

type ODETest
    ivp::IVPWithSolution
    solver::Function
    rtol::Real
    atol::Real
end

# Single equations

# A1: The negative potential
A1 = IVPWithSolution( "A1", (t,y)->-y , 1. , [0., 10.] , t->exp(-t) )
# A2: A special case of the Riccati equation
A2 = IVPWithSolution( "A1", (t,y)->-y^3 / 2, 1, [0, 5], t-> 1/sqrt(t+1) )



end