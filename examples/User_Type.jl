#=
This is an example of how to use a custom user type of with ODE module.

An additional example can be found in `test/inverface-tests.jl`

The ODE considered in this example is the simple oscilator

  a' = -b
  b' =  a

with initial condtion a(t=0) = 0 and b(t=0) = 0.1.
=#

using ODE
using PyPlot
using Compat

const DataFloat = Float64
mutable struct Vec2
    a::DataFloat
    b::DataFloat
end

# This is the minimal set of function required on the type inorder to work with
# the ODE module
Base.:/(x::Vec2, y::Real) = Vec2(x.a / y, x.b / y)
Base.:*(y::Real, x::Vec2) = Vec2(y * x.a, y * x.b)
Base.:*(x::Vec2, y::Real) = y * x
Base.:.*(y::Real, x::Vec2) = y * x
Base.:+(x::Vec2, y::Vec2) = Vec2(x.a + y.a, x.b + y.b)
Base.:-(x::Vec2, y::Vec2) = Vec2(x.a - y.a, x.b - y.b)
Base.norm(x::Vec2) = sqrt(x.a^2 + x.b^2)
Base.zero(x::Type{Vec2}) = Vec2(zero(DataFloat), zero(DataFloat))
ODE.isoutofdomain(x::Vec2) = isnan(x.a) || isnan(x.b)

# RHS function
f(t, y) = Vec2(-y.b, y.a)

# Initial condtions
start = Vec2(0.0, 0.1)

# Time vector going from 0 to 2*PI in 0.01 increments
time = 0:0.1:(4 * pi)

# Solve the ODE
t, y = ode45(f, start, time)

# Plot the solution
a = map(y -> y.a, y)
b = map(y -> y.b, y)
plot(t, a, label = "a(t)")
plot(t, b, label = "b(t)")
legend()
