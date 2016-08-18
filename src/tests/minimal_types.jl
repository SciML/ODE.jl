import Base: <, >, >=, <=
import Base: +, -, *, /, ^
import Base: .+, .-
import Base: copy, zero, getindex, setindex!, similar
import Base: eps, convert, promote_rule


# position variable
type Position{T} <: AbstractVector{T}
    x::T
    y::T
end

copy(p::Position) = Position(p.x,p.y)
similar{T}(p::Position{T}) = copy(p)

for op in (:+, :-, :.+)
    @eval ($op)(p1::Position,p2::Position) = Position(($op)(p1.x,p2.x),
                                                      ($op)(p1.y,p2.y))
end

Base.size(::Position)=(2,)
getindex{T}(p::Position{T},i::Int)         = i==1 ? p.x : p.y
setindex!{T}(p::Position{T},val::T,i::Int) = i==1 ? p.x=val : p.y=val

# MyFloat variable
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

# binary operators
for op in (:+, :-, :*, :/, :^)
    @eval ($op)(t1::MyFloat,t2::MyFloat) = MyFloat(($op)(t1.t,t2.t))
end

# unary operators
for op in (:-,)
    @eval ($op)(t::MyFloat) = MyFloat(($op)(t.t))
end

# comparisons
for op in (:<, :>, :>=, :<=)
    @eval ($op)(t1::MyFloat,t2::MyFloat) = ($op)(t1.t,t2.t)
end

# vector times scalar multiplication
*(t::MyFloat,p::Position)=Position(t*p.x,t*p.y)
*(p::Position,t::MyFloat)= *(t,p)


const case_minimal_type =
    Dict(:ivp      => ExplicitODE(MyFloat(0.0),
                                  Position(MyFloat(0.0),MyFloat(1.0)),
                                  (t,y,dy)->(dy.x=y.y;dy.y=-y.x)),
         :sol      => t->Position(sin(t),cos(t)),
         :name     => "y'=6 (minimal types)",
         :options  => Dict(:initstep => MyFloat(0.1),
                           :tstop => MyFloat(1.0))
         )

Base.sin{T}(t::MyFloat{T})=MyFloat{T}(sin(t.t))
Base.cos{T}(t::MyFloat{T})=MyFloat{T}(cos(t.t))
