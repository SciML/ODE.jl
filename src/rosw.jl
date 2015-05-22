# Rosenbrock-Wanner methods
###########################

export ode_rosw, ode_rosw_fixed

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
ode_rosw_fixed(fn, Jfn, y0, tspan) = oderosw_fixed(fn, Jfn, y0, tspan, bt_ros34pw2)
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

ode_rosw(fn, Jfn, y0, tspan;kwargs...) = oderosw_adapt(fn, Jfn, y0, tspan, bt_ros34pw2; kwargs...)
function oderosw_adapt{N,S}(fn, Jfn, y0::AbstractVector, tspan, btab::TableauRosW{N,S};
                            reltol = 1.0e-5, abstol = 1.0e-8,
                            norm=Base.norm,
                            minstep=abs(tspan[end] - tspan[1])/1e9,
                            maxstep=abs(tspan[end] - tspan[1])/2.5,
                            initstep=0.,
                            points=:all
                            )
    fn_expl = (t,y)->fn(y, y*0)
    Et, Eyf, Ty, btab = make_consistent_types(fn, y0, tspan, btab)
    btab = tabletransform(btab)
    # parameters
    order = minimum(btab.order)
    timeout_const = 5 # after step reduction do not increase step for
                      # timeout_const steps

    ## Initialization
    dof = length(y0)
    tspan = convert(Vector{Et}, tspan)
    t = tspan[1]
    tstart = tspan[1]
    tend = tspan[end]

    # work arrays:
    y      = similar(y0, Eyf, dof)      # y at time t
    y[:]   = y0
    ytrial = similar(y0, Eyf, dof) # trial solution at time t+dt
    yerr   = similar(y0, Eyf, dof) # error of trial solution
    ks = Array(Ty, S)
    # allocate!(ks, y0, dof) # no need to allocate as fn is not in-place
    ytmp   = similar(y0, Eyf, dof)

    # output ys
    nsteps_fixed = length(tspan) # these are always output
    ys = Array(Ty, nsteps_fixed)
    allocate!(ys, y0, dof)
    ys[1] = y0

    # Option points determines where solution is returned:
    if points==:all
        tspan_fixed = tspan
        tspan = Et[tstart]
        iter_fixed = 2 # index into tspan_fixed
        sizehint!(tspan, nsteps_fixed)
    elseif points!=:specified
        error("Unrecognized option points==$points")
    end
    # Time
    dt, tdir, ks[1] = hinit(fn_expl, y, tstart, tend, order, reltol, abstol) # sets ks[1]=f0
    if initstep!=0
        dt = sign(initstep)==tdir ? initstep : error("initstep has wrong sign.")
    end
    # Diagnostics
    dts = Et[]
    errs = Float64[]
    steps = [0,0]  # [accepted, rejected]

    ## Integration loop
    laststep = false
    timeout = 0 # for step-control
    iter = 2 # the index into tspan and ys
    while true
        ytrial[:] = rosw_step(fn, Jfn, y, dt, btab, 1)
        yerr[:] = ytrial - rosw_step(fn, Jfn, y, dt, btab, 2)
        err, newdt, timeout = stepsize_hw92!(dt, tdir, y, ytrial, yerr, order, timeout,
                                            dof, abstol, reltol, maxstep, norm)

        if err<=1.0 # accept step
            # diagnostics
            steps[1] +=1
            push!(dts, dt)
            push!(errs, err)

            # Output:
            f0 = ks[1]
            f1 = fn_expl(t+dt, ytrial)
            if points==:specified
                # interpolate onto given output points
                while iter-1<nsteps_fixed && (tdir*tspan[iter]<tdir*(t+dt) || laststep) # output at all new times which are < t+dt
                    hermite_interp!(ys[iter], tspan[iter], t, dt, y, ytrial, f0, f1) # TODO: 3rd order only!
                    iter += 1
                end
            else
                # first interpolate onto given output points
                while iter_fixed-1<nsteps_fixed && tdir*t<tdir*tspan_fixed[iter_fixed]<tdir*(t+dt) # output at all new times which are < t+dt
                    yout = hermite_interp(tspan_fixed[iter_fixed], t, dt, y, ytrial, f0, f1)
                    index_or_push!(ys, iter, yout) # TODO: 3rd order only!
                    push!(tspan, tspan_fixed[iter_fixed])
                    iter_fixed += 1
                    iter += 1
                end
                # but also output every step taken
                index_or_push!(ys, iter, copy(ytrial))
                push!(tspan, t+dt)
                iter += 1
            end
            ks[1] = f1 # load ks[1]==f0 for next step

            # Break if this was the last step:
            laststep && break

            # Swap bindings of y and ytrial, avoids one copy
            y, ytrial = ytrial, y

            # Update t to the time at the end of current step:
            t += dt
            dt = newdt

            # Hit end point exactly if next step within 1% of end:
            if tdir*(t+dt*1.01) >= tdir*tend
                dt = tend-t
                laststep = true # next step is the last, if it succeeds
            end
        elseif abs(newdt)<minstep  # minimum step size reached, break
            println("Warning: dt < minstep.  Stopping.")
            break
        else # redo step with smaller dt
            laststep = false
            steps[2] +=1
            dt = newdt
            timeout = timeout_const
        end
    end
    return tspan, ys

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
