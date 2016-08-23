import Base: <, >, >=, <=
import Base: +, -, *, /, ^
# for the vector type
import Base: getindex, setindex!, similar
# for the scalar type
import Base: eps, convert, promote_rule, sqrt, isfinite


# position variable
type Position{T} <: AbstractVector{T}
    x::T
    y::T
end

similar{T}(p::Position{T},::Type{T}) = Position(p.x,p.y)

for op in (:+, :-)
    @eval ($op)(p1::Position,p2::Position) = Position(($op)(p1.x,p2.x),
                                                      ($op)(p1.y,p2.y))
end

Base.size(::Position)=(2,)
getindex{T}(p::Position{T},i::Int)         = i==1 ? p.x : p.y
setindex!{T}(p::Position{T},val::T,i::Int) = i==1 ? p.x=val : p.y=val

# MyFloat variable.  It can be constructed from any type but cannot be
# converted, so operations between MyFloat{Float64} and Float64 should
# throw a conversion error.  This is to detect operations mixing high
# precision floats (like BigFloat) with lower precision constants,
# which could result in decreasing the overall precision of the
# algorithm.
immutable MyFloat <: Real
    t::Float64
    # the dummy is here to prevent the use of the default constructor
    dummy::Int
end

# we need these to construct MyFloat from constants, constants are
# predefined in terms of numbers of inifinite precision such as Int or
# Rational
convert(::Type{MyFloat},s::Rational) = MyFloat(convert(Float64,s),0)
convert(::Type{MyFloat},s::Integer) = MyFloat(convert(Float64,s),0)
promote_rule{R<:Rational}(::Type{MyFloat},::Type{R}) = MyFloat
promote_rule{R<:Integer}(::Type{MyFloat},::Type{R}) = MyFloat

eps(::Type{MyFloat}) = MyFloat(eps(Float64),0)
isfinite(x::MyFloat) = isfinite(x.t)
# necessary for the modified Rosenbrock integrator
sqrt(t::MyFloat)=MyFloat(sqrt(t.t),0)

# binary operators
for op in (:+, :-, :*, :/, :^)
    @eval ($op)(t1::MyFloat,t2::MyFloat) = MyFloat(($op)(t1.t,t2.t),0)
end

# See #18114
^(x::MyFloat, y::Rational) = x^(convert(MyFloat,y.num)/convert(MyFloat,y.den))

# unary operators
for op in (:-,)
    @eval ($op)(t::MyFloat) = MyFloat(($op)(t.t),0)
end

# comparisons
for op in (:<, :>, :>=, :<=)
    @eval ($op)(t1::MyFloat,t2::MyFloat) = ($op)(t1.t,t2.t)
end

# these are only necessary because they are used in the definition of
# the ODE test case (see test_cases.jl, :harmonic_minimal_types)
Base.sin(t::MyFloat)=MyFloat(sin(t.t),0)
Base.cos(t::MyFloat)=MyFloat(cos(t.t),0)
