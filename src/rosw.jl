# Rosenbrock-Wanner methods
###########################

export ode_rosw

# Rosenbrock-W methods are typically specified for autonomous DAE:
# Mẋ = f(x)
#
# by the stage equations
# M kᵢ = hf(x₀ + + Σⱼ aᵢⱼkⱼ) + h J Σⱼ γᵢⱼkⱼ
#
# and step completion formula
# x₁ = x₀ + \Sigmaⱼ bⱼ kⱼ
#
# The method used here uses transformed equations as done in PETSc.


# Tableaus
##########

immutable TableauRosW{Name, S, T} <: Tableau{Name, S, T}
    order::(@compat(Tuple{Vararg{Int}})) # the order of the methods
    a::Matrix{T}
    γ::Matrix{T}
    # one or several row vectors.  First row is used for the step,
    # second for error calc.
    b::Matrix{T}
    function TableauRosW(order,a,γ,b)
        @assert isa(S,Integer)
        @assert isa(Name,Symbol)
        @assert S==size(γ,1)==size(a,2)==size(γ,1)==size(a,2)==size(b,2)
        @assert size(b,1)==length(order)
        new(order,a,γ,b)
    end
end
function TableauRosW{T}(name::Symbol, order::(@compat(Tuple{Vararg{Int}})),
                   a::Matrix{T}, γ::Matrix{T}, b::Matrix{T})
    TableauRosW{name,size(b,2),T}(order, a, γ, b)
end
function TableauRosW(name::Symbol, order::(@compat(Tuple{Vararg{Int}})), T::Type,
                   a::Matrix, γ::Matrix, b::Matrix)
    TableauRosW{name,size(b,2),T}(order, convert(Matrix{T},a),
                                  convert(Matrix{T},γ), convert(Matrix{T},b) )
end
conv_field{T,N}(D,a::Array{T,N}) = convert(Array{D,N}, a)
function Base.convert{Tnew<:Real,Name,S,T}(::Type{Tnew}, tab::TableauRosW{Name,S,T})
    # Converts the tableau coefficients to the new type Tnew
    newflds = ()
    @compat for n in fieldnames(tab)
        fld = getfield(tab,n)
        if eltype(fld)==T
            newflds = tuple(newflds..., conv_field(Tnew, fld))
        else
            newflds = tuple(newflds..., fld)
        end
    end
    TableauRosW{Name,S,Tnew}(newflds...) # TODO: could this be done more generically in a type-stable way?
end



# Transformed Tableau, used only internally
immutable TableauRosW_T{Name, S, T} <: Tableau{Name, S, T}
    order::(@compat(Tuple{Vararg{Int}})) # the order of the methods
    a::Matrix{T}  # this is TableauRosW.a transformed
    γii::Vector{T}
    γinv::Matrix{T}
    b::Matrix{T} # this is TableauRosW.b transformed
end
function tabletransform{Name,S,T}(rt::TableauRosW{Name,S,T})
    γii = diag(rt.γ)
    γinv = inv(rt.γ)
    ahat = rt.a * γinv
    bhat = similar(rt.b)
    bhat[1,:] = squeeze(rt.b[1,:]*γinv,1)
    bhat[2,:] = squeeze(rt.b[2,:]*γinv,1)
    return TableauRosW_T{Name,S,T}(rt.order, ahat, γii, γinv, bhat)
end


## tableau for ros34pw2
const bt_ros34pw2 = TableauRosW(:ros34pw2, (3,4), Float64,
                                [0  0  0  0
                                 8.7173304301691801e-01  0  0  0
                                 8.4457060015369423e-01  -1.1299064236484185e-01  0  0
                                 0  0  1.  0],
                                [4.3586652150845900e-01  0  0  0
                                 -8.7173304301691801e-01  4.3586652150845900e-01  0  0
                                 -9.0338057013044082e-01  5.4180672388095326e-02  4.3586652150845900e-01  0
                                 2.4212380706095346e-01  -1.2232505839045147e+00  5.4526025533510214e-01  4.3586652150845900e-01],
                                 [2.4212380706095346e-01  -1.2232505839045147e+00  1.5452602553351020e+00  4.3586652150845900e-01
                                  3.7810903145819369e-01  -9.6042292212423178e-02  5.0000000000000000e-01  2.1793326075422950e-01]
                                  )

###################
# Fixed step solver
###################
ode_rosw(fn, Jfn, y0, tspan) = oderosw_fixed(fn, Jfn, y0::AbstractVector, tspan, bt_ros34pw2)
function oderosw_fixed{N,S}(fn, Jfn, y0::AbstractVector, tspan,
                            btab::TableauRosW{N,S})
    # TODO: refactor with oderk_fixed
    Et, Eyf, Ty, btab = make_consistent_types(fn, y0, tspan, btab)
    btab = tabletransform(btab)
    dof = length(y0)

    ys = Array(Ty, length(tspan))
    allocate!(ys, y0, dof)
    ys[1] = deepcopy(y0)

    tspan = convert(Vector{Et}, tspan)
    # work arrays:
    ks = Array(Ty, S)
    # allocate!(ks, y0, dof) # no need to allocate as fn is not in-place
    ytmp = similar(y0, Eyf, dof)
    for i=1:length(tspan)-1
        dt = tspan[i+1]-tspan[i]
        ys[i+1] = rosw_step(fn, Jfn, ys[i], dt, btab, 2)
    end
    return tspan, ys
end

function linsolve(h, hprime, y0)
    # Does one Newton step, i.e. a linear solve.
    return y0 - hprime(y0)\h(y0)
end


function oderosw_adapt()

end

function rosw_step(g, gprime, x0, dt,
                   btab::TableauRosW_T, bt_ind=1)
    # This takes one step for a ode/dae system defined by
    # g(x,xdot)=0
    # gprime(x, xdot, α) = dg/dx + α dg/dxdot

    stages = size(btab.a,1)
    dof = length(x0)
    ys  = zeros(eltype(x0), stages, dof)

    # stage solutions
    jacobian_stale = true
    hprime_store = zeros(dof,dof)

    cij = (tril(btab.γinv)-diagm(diag(btab.γinv)))/dt
    for i=1:stages
        u = zeros(dof)
        udot = zeros(dof)
        function h(yi) # length(yi)==dof
            # this function is Eq 5
            u[:] = x0 + yi
            for j=1:i-1
                for d=1:dof
                    u[d] += btab.a[i,j]*ys[j,d]
                end
            end
            udot[:] = 1./(dt*btab.γii[i]).*yi  # jed: is this index with the stage?
            for j=1:i
                # this is Eq.5-3
                for d=1:dof
                    udot[d] += cij[i,j]*ys[j,d]
                end
            end
            g(u, udot)
        end
        function hprime(yi)
            # here we only update the jacobian once per time step
            if jacobian_stale
                hprime_store[:,:] = gprime(u, udot, 1./(dt*btab.γii[i]))  # jed: is this index with the stage?
            end
            hprime_store
        end
        # this is just a linear solve, usually
#        ys[i,:] = nlsolve(h, zeros(dof), hprime; opts=one_step_only)
        ys[i,:] = linsolve(h, hprime, zeros(dof))
#        ys[i,:] = newtonsolve(h, hprime, zeros(dof), maxsteps=1, warn=false)
        if i==1
            # calculate jacobian
            jacobian_stale = false
        end
    end

    # completion:  (done twice for error control)
    x1 = zeros(eltype(x0), length(x0))
    for i=1:dof
        for j=1:stages
            x1[i] += btab.b[bt_ind,j]*ys[j,i]
        end
    end
    return x0 + x1
end
