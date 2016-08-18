import Base: <, >, >=, <=
import Base: +, -, *, /, ^
# for the vector type
import Base: getindex, setindex!, similar
# for the scalar type
import Base: eps, convert, promote_rule, sqrt


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
immutable MyFloat{T} <: Real
    t::T
end

# we need these to construct MyFloat from constants, constants are
# predefined in terms of numbers of inifinite precision such as Int or
# Rational
convert{T}(::Type{MyFloat{T}},s::Rational) = MyFloat{T}(convert(T,s))
convert{T}(::Type{MyFloat{T}},s::Integer) = MyFloat{T}(convert(T,s))
promote_rule{T<:MyFloat,R<:Rational}(::Type{T},::Type{R}) = T
promote_rule{T<:MyFloat,R<:Integer}(::Type{T},::Type{R}) = T

eps{T}(::Type{MyFloat{T}}) = MyFloat{T}(eps(T))
# necessary for the modified Rosenbrock integrator
sqrt{T}(t::MyFloat{T})=MyFloat{T}(sqrt(t.t))

# binary operators
for op in (:+, :-, :*, :/, :^)
    @eval ($op){T}(t1::MyFloat{T},t2::MyFloat{T}) = MyFloat{T}(($op)(t1.t,t2.t))
end

# See #18114
^{T<:MyFloat}(x::T, y::Rational) = x^(convert(T,y.num)/convert(T,y.den))

# unary operators
for op in (:-,)
    @eval ($op)(t::MyFloat) = MyFloat(($op)(t.t))
end

# comparisons
for op in (:<, :>, :>=, :<=)
    @eval ($op){T}(t1::MyFloat{T},t2::MyFloat{T}) = ($op)(t1.t,t2.t)
end

# these are only necessary because they are used in the definition of
# the ODE test case (see test_cases.jl, :harmonic_minimal_types)
Base.sin{T}(t::MyFloat{T})=MyFloat{T}(sin(t.t))
Base.cos{T}(t::MyFloat{T})=MyFloat{T}(cos(t.t))
