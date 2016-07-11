using ODE
using PyPlot

typealias DataFloat Float64

type Vec2
  a::DataFloat
  b::DataFloat
end

# Necessary for the solver
Base.(:/)(x::Vec2, y::Real) = Vec2(x.a/y, x.b/y)
Base.(:*)(y::Real, x::Vec2) = Vec2(y * x.a, y * x.b)
Base.(:*)(x::Vec2, y::Real) = y*x
Base.(:.*)(y::Real, x::Vec2) = y*x
Base.(:+)(x::Vec2, y::Vec2) = Vec2(x.a + y.a, x.b + y.b)
Base.(:-)(x::Vec2, y::Vec2) = Vec2(x.a - y.a, x.b - y.b)
Base.norm(x::Vec2) = sqrt(x.a^2 + x.b^2)
Base.zero(x::Type{Vec2}) = Vec2(zero(DataFloat), zero(DataFloat))
Base.isnan(x::Vec2) = isnan(x.a) || isnan(x.b)

function f(t, y)
  (a, b) = (y.a, y.b)
  Vec2(-b, a)
end

let
  # Initial condtions
  const start = Vec2(0.0, 0.1)

  # Time vector going from 0 to 2*PI in 0.01 increments
  const time = 0:0.1:4*pi

  # Solve the ODE
  t, y = ode45(f, start, time)

  # Plot the solution
  a = map(y -> y.a, y)
  b = map(y -> y.b, y)
  plot(t, a, label="x(t)")
  plot(t, b, label="x'(t)")
  legend()
end
