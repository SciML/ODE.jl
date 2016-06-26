###############################################################################
## Coefficient Tableaus
###############################################################################

# Butcher Tableaus, or more generally coefficient tables
# see Hairer & Wanner 1992, p. 134, 166

abstract Tableau{T<:Real}
# Name is the name of the tableau/method (a symbol)
# S is the number of stages (an int)
# T is the type of the coefficients
#
# TODO: have a type parameter which specifies adaptive vs non-adaptive
#
# For all types of tableaus it assumes fields:
# order::(Int...) # order of the method(s)
#
# For Runge-Kutta methods it assumes fields:
# a::Matrix{T}  # SxS matrix
# b::Matrix{T}  # 1 or 2 x S matrix (fixed step/ adaptive)
# c::Vector{T}  # S
#
# For a tableau:
#  c1  | a_11   ....   a_1s
#  .   | a_21 .          .
#  .   | a_31     .      .
#  .   | ....         .  .
#  c_s | a_s1  ....... a_ss
# -----+--------------------
#      | b_1     ...   b_s   this is the one used for stepping
#      | b'_1    ...   b'_s  this is the one used for error-checking

Base.eltype{T}(b::Tableau{T}) = T
order(b::Tableau) = b.order
# Subtypes need to define a convert method to convert to a different
# eltype with signature:
Base.convert{Tnew<:Real}(::Type{Tnew}, tab::Tableau) = error("Define convert method for concrete Tableau types")


###########################################
# Tableaus for explicit Runge-Kutta methods
###########################################

# m3: these tableaus should go into rk.jl as they belong to R-K methods

immutable TableauRKExplicit{T} <: Tableau{T}
    order::(@compat(Tuple{Vararg{Int}})) # the order of the methods
    a::Matrix{T}
    # one or several row vectors.  First row is used for the step,
    # second for error calc.
    b::Matrix{T}
    c::Vector{T}
    isFSAL::Bool
    s::Int
    name :: AbstractString
    function TableauRKExplicit(name,order,a,b,c)
        s = length(c)
        @assert c[1]==0
        @assert istril(a)
        @assert s==size(a,1)==size(a,2)==size(b,2)
        @assert size(b,1)==length(order)
        @assert norm(sum(a,2)-c'',Inf)<T(1e-10) # consistency.
        isFSAL = (a[end,:]==b[1,:] && c[end]==1)
        new(order,a,b,c,isFSAL,s,name)
    end
end


function TableauRKExplicit{T}(name::AbstractString, order::(@compat(Tuple{Vararg{Int}})),
                              a::Matrix{T}, b::Matrix{T}, c::Vector{T})
    TableauRKExplicit{T}(name, order, a, b, c)
end


function TableauRKExplicit(name::AbstractString, order::(@compat(Tuple{Vararg{Int}})), T::Type,
                           a::Matrix, b::Matrix, c::Vector)
    TableauRKExplicit{T}(name, order, convert(Matrix{T},a),
                         convert(Matrix{T},b), convert(Vector{T},c) )
end


# TODO: remove conv_field
conv_field{T,N}(D,a::Array{T,N}) = convert(Array{D,N}, a)


lengthks(tab::TableauRKExplicit) = length(tab.c)

# TODO: there should be a better way to do it
function Base.convert{Tnew<:Real,T}(::Type{Tnew}, tab::TableauRKExplicit{T})
    # Converts the tableau coefficients to the new type Tnew
    newflds = ()
    @compat for n in [:a,:b,:c]
        fld = getfield(tab,n)
        if eltype(fld)==T
            newflds = tuple(newflds..., conv_field(Tnew, fld))
        else
            newflds = tuple(newflds..., fld)
        end
    end
    TableauRKExplicit{Tnew}(tab.name, tab.order, newflds...) # TODO: could this be done more generically in a type-stable way?
end


isexplicit(b::TableauRKExplicit) = istril(b.a) # Test whether it's an explicit method
isadaptive(b::TableauRKExplicit) = size(b.b, 1)==2

## Tableaus for explicit RK methods
# Fixed step:
const tableaus_rk_explicit = Dict{Symbol,ODE.TableauRKExplicit{Rational{Int}}}()

tableaus_rk_explicit[:feuler] =
    TableauRKExplicit("Forward Euler",(1,), Rational{Int64},
                      zeros(Int,1,1),
                      [1]',
                      [0]
                      )

tableaus_rk_explicit[:midpoint] =
    TableauRKExplicit("Midpoint",(2,), Rational{Int64},
                      [0  0
                       1//2  0],
                      [0, 1]',
                      [0, 1//2]
                      )

tableaus_rk_explicit[:heun] =
    TableauRKExplicit("Heun",(2,), Rational{Int64},
                      [0  0
                       1  0],
                      [1//2, 1//2]',
                      [0, 1])

tableaus_rk_explicit[:rk4] =
    TableauRKExplicit("Runge-Kutta(4)",(4,),Rational{Int64},
                      [0    0    0 0
                       1//2 0    0 0
                       0    1//2 0 0
                       0    0    1 0],
                      [1//6, 1//3, 1//3, 1//6]',
                      [0, 1//2, 1//2, 1])

# Adaptive step:
# Heun Euler https://en.wikipedia.org/wiki/Runge–Kutta_methods
tableaus_rk_explicit[:rk21] =
    TableauRKExplicit("Heun Euler",(2,1), Rational{Int64},
                      [0     0
                       1     0],
                      [1//2  1//2
                       1     0],
                      [0,    1])

# Bogacki–Shampine coefficients
tableaus_rk_explicit[:rk23] =
    TableauRKExplicit("Bogacki-Shampine",(2,3), Rational{Int64},
                      [0           0      0      0
                       1/2         0      0      0
                       0         3/4      0      0
                       2/9       1/3     4/9     0],
                      [7/24 1/4 1/3 1/8
                       2/9 1/3 4/9 0],
                      [0, 1//2, 3//4, 1]
                      )

# Fehlberg https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method
tableaus_rk_explicit[:rk45] =
    TableauRKExplicit("Fehlberg",(4,5),Rational{Int64},
                      [  0           0           0            0         0     0
                         1//4        0           0            0         0     0
                         3//32       9//32       0            0         0     0
                         1932//2197 -7200//2197  7296//2197      0         0     0
                         439//216     -8        3680//513    -845//4104   0     0
                         -8//27       2       -3544//2565   1859//4104 -11//40 0 ],
[25//216      0        1408//2565   2197//4104  -1//5  0
 16//135      0        6656//12825 28561//56430 -9//50 2//55],
[0,          1//4,       3//8,       12//13,    1,    1//2])

# Dormand-Prince https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method
tableaus_rk_explicit[:dopri5] =
    TableauRKExplicit("Dormand-Prince", (5,4), Rational{Int64},
                     [0           0            0                   0            0      0 0
                      1//5        0            0                   0            0      0 0
                      3//40       9//40        0                   0            0      0 0
                      44//45      -56//15      32//9               0            0      0 0
                      19372//6561 -25360//2187 64448//6561 -212//729            0      0 0
                      9017//3168  -355//33     46732//5247   49//176 -5103//18656      0 0
                      35//384     0            500//1113    125//192  -2187//6784 11//84 0],
                     [35//384         0     500//1113       125//192      -2187//6784         11//84      0
                      5179//57600     0     7571//16695     393//640     -92097//339200     187//2100     1//40],
                     [0, 1//5, 3//10, 4//5, 8//9, 1, 1]
                     )

# Fehlberg 7(8) coefficients
# Values from pag. 65, Fehlberg, Erwin. "Classical fifth-, sixth-, seventh-, and eighth-order Runge-Kutta formulas with stepsize control".
# National Aeronautics and Space Administration.
tableaus_rk_explicit[:feh78] =
    TableauRKExplicit("Fehlberg(7,8)", (7,8), Rational{Int64},
                            [     0       0      0       0         0          0        0        0      0       0     0 0 0
                                  2//27   0      0       0         0          0        0        0      0       0     0 0 0
                                  1//36   1//12  0       0         0          0        0        0      0       0     0 0 0
                                  1//24   0      1//8    0         0          0        0        0      0       0     0 0 0
                                  5//12   0    -25//16  25//16     0          0        0        0      0       0     0 0 0
                                  1//20   0      0       1//4      1//5       0        0        0      0       0     0 0 0
                                -25//108  0      0     125//108  -65//27    125//54    0        0      0       0     0 0 0
                                 31//300  0      0       0       61//225    -2//9     13//900   0      0       0     0 0 0
                                  2       0      0     -53//6    704//45   -107//9    67//90    3      0       0     0 0 0
                                -91//108  0      0      23//108 -976//135   311//54  -19//60   17//6  -1//12   0     0 0 0
                               2383//4100 0      0    -341//164 4496//1025 -301//82 2133//4100 45//82 45//164 18//41 0 0 0
                                  3//205  0      0       0        0        -6//41     -3//205  -3//41  3//41   6//41 0 0 0
                              -1777//4100 0      0    -341//164 4496//1025 -289//82 2193//4100 51//82 33//164 12//41 0 1 0],
                              [41//840 0 0 0 0 34//105 9//35 9//35 9//280 9//280 41//840 0 0
                               0 0 0 0 0 34//105 9//35 9//35 9//280 9//280 0     41//840 41//840],
                               [0,    2//27, 1//9, 1//6 , 5//12, 1//2 , 5//6 , 1//6 , 2//3 , 1//3 , 1 , 0, 1]
                            )
