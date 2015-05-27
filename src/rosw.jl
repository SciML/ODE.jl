# Rosenbrock-Wanner methods
###########################
#
# Main references:
# - Wanner & Hairer 1996
# - PETSc http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSROSW.html
export ode_rosw, ode_rosw_fixed

# TODO:
# - AD Jacobian
# - fix fixed step solver

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
    γinv::Matrix{T}
    b::Matrix{T} # this is TableauRosW.b transformed
    # derived quantities:
    γii::T
    c::Matrix{T} # = (tril(btab.γinv)-diagm(diag(btab.γinv)))
end
function tabletransform{Name,S,T}(rt::TableauRosW{Name,S,T})
    # the code only works if
    if !all(x->x==rt.γ[1],diag(rt.γ))
        error("This Rosenbrock implementation only works for tableaus with γ_ii==γ_jj for all i,j.")
    end
    γii = rt.γ[1,1]
    γinv = inv(rt.γ)
    ahat = rt.a * γinv
    bhat = similar(rt.b)
    bhat[1,:] = squeeze(rt.b[1,:]*γinv,1)
    bhat[2,:] = squeeze(rt.b[2,:]*γinv,1)
    c = (tril(γinv)-diagm(diag(γinv))) # negative of W&H definition
    return TableauRosW_T{Name,S,T}(rt.order, ahat, γinv, bhat, γii, c)
end


## tableau for ros34pw2 http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSROSWRA34PW2.html
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

## tableau for RODAS3 http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/TS/TSROSWRODAS3.html
const bt_ros_rodas3 = TableauRosW(:ros_rodas3, (3,4), Rational{Int},
                                  [0 0 0 0
                                   0 0 0 0
                                   1 0 0 0
                                   3//4 -1//4 1//2 0],
                                  [1//2 0 0 0
                                   1  1//2 0 0
                                   -1//4 -1//4 1//2 0
                                   1//12 1//12 -2//3 1//2],
                                   [5//6 -1//6 -1//6 1//2
                                    3//4 -1//4 1//2 0]
                                    )

###################
# Fixed step solver
###################
ode_rosw_fixed(fn, x0, tspan; kwargs...) = oderosw_fixed(fn, x0, tspan, bt_ros34pw2, kwargs...)
function oderosw_fixed{N,S}(fn, x0::AbstractVector, tspan,
                            btab::TableauRosW{N,S};
                            jacobian=numerical_jacobian(fn, 1e-6, 1e-6)
                            )
    # TODO: refactor with oderk_fixed
    Et, Exf, Tx, btab = make_consistent_types(fn, x0, tspan, btab)
    btab = tabletransform(btab)
    dof = length(x0)

    xs = Array(Tx, length(tspan))
    allocate!(xs, x0, dof)
    xs[1] = deepcopy(x0)

    tspan = convert(Vector{Et}, tspan)
    # work arrays:
    k = Array(Tx, S) # stage variables
    allocate!(k, x0, dof)
    ks = zeros(Exf,dof) # work vector for one k
    jac_store = zeros(Exf,dof,dof) # Jacobian storage
    u = zeros(Exf,dof)    # work vector 
    udot = zeros(Exf,dof) # work vector 
    
    # allocate!(ks, x0, dof) # no need to allocate as fn is not in-place
    xtmp = similar(x0, Exf, dof)
    for i=1:length(tspan)-1
        dt = tspan[i+1]-tspan[i]
        rosw_step!(xs[i+1], fn, jacobian, xs[i], dt, dof, btab,
                   k, jac_store, ks, u, udot)
    end
    return tspan, xs
end


ode_rosw(fn, x0, tspan;kwargs...) = oderosw_adapt(fn, x0, tspan, bt_ros34pw2; kwargs...)
ode_rodas3(fn, x0, tspan;kwargs...) = oderosw_adapt(fn, x0, tspan, bt_ros_rodas3; kwargs...)
function oderosw_adapt{N,S}(fn, x0::AbstractVector, tspan, btab::TableauRosW{N,S};
                            reltol = 1.0e-5, abstol = 1.0e-8,
                            norm=Base.norm,
                            minstep=abs(tspan[end] - tspan[1])/1e18,
                            maxstep=abs(tspan[end] - tspan[1])/2.5,
                            initstep=0.,
                            jacobian=numerical_jacobian(fn, reltol, abstol)
#                            points=:all
                            )
    # TODO: refactor with oderk_adapt

    # FIXME: add interpolation
    if length(tspan)>2
        error("specified output times not supported yet")
    else
        points=:all
    end
    ## Figure types
    fn_expl = (t,x)->(out=similar(x); fn(out, x, x*0); out)
    Et, Exf, Tx, btab = make_consistent_types(fn, x0, tspan, btab)
    btab = tabletransform(btab)
    # parameters
    order = minimum(btab.order)
    timeout_const = 5 # after step reduction do not increase step for
                      # timeout_const steps

    ## Setup
    dof = length(x0)
    tspan = convert(Vector{Et}, tspan)
    t = tspan[1]
    tstart = tspan[1]
    tend = tspan[end]

    # work arrays:
    x      = similar(x0, Exf, dof) # x at time t (time at beginning of step)
    x[:]   = x0                    # fill with IC
    xtrial = similar(x0, Exf, dof) # trial solution at time t+dt
    xerr   = similar(x0, Exf, dof) # error of trial solution
    k = Array(Tx, S) # stage variables
    allocate!(k, x0, dof)
    ks = zeros(Exf,dof) # work vector for one k
    jac_store = zeros(Exf,dof,dof) # Jacobian storage
    u = zeros(Exf,dof)    # work vector 
    udot = zeros(Exf,dof) # work vector 

    # output xs
    nsteps_fixed = length(tspan) # these are always output
    xs = Array(Tx, nsteps_fixed)
    allocate!(xs, x0, dof)
    xs[1] = x0

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
    dt, tdir, f0 = hinit(fn_expl, x, tstart, tend, order, reltol, abstol)
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
    iter = 2 # the index into tspan and xs
    while true
        rosw_step!(xtrial, fn, jacobian, x, dt, dof, btab,
                   k, jac_store, ks, u, udot)
        # Completion again for embedded method, see line 927 of
        # http://www.mcs.anl.gov/petsc/petsc-current/src/ts/impls/rosw/rosw.c.html#TSROSW
        xerr[:] = 0.0
        for s=1:S
            for d=1:dof
                xerr[d] += (btab.b[2,s]-btab.b[1,s])*k[s][d]
            end
        end
        err, newdt, timeout = stepsize_hw92!(dt, tdir, x, xtrial, xerr, order, timeout,
                                            dof, abstol, reltol, maxstep, norm)
        if err<=1.0 # accept step
            # diagnostics
            steps[1] +=1
            push!(dts, dt)
            push!(errs, err)

            # Output:
            # f0 = k[1]
            # f1 = fn_expl(t+dt, xtrial)
            if points==:specified
                # interpolate onto given output points
                while iter-1<nsteps_fixed && (tdir*tspan[iter]<tdir*(t+dt) || laststep) # output at all new times which are < t+dt
                    error("FIXME: not implemented yet")
                    hermite_interp!(xs[iter], tspan[iter], t, dt, x, xtrial, f0, f1) # TODO: 3rd order only!
                    iter += 1
                end
            else
                # first interpolate onto given output points
                while iter_fixed-1<nsteps_fixed && tdir*t<tdir*tspan_fixed[iter_fixed]<tdir*(t+dt) # output at all new times which are < t+dt
                    error("FIXME: not implemented yet")
                    xout = hermite_interp(tspan_fixed[iter_fixed], t, dt, x, xtrial, f0, f1)
                    index_or_push!(xs, iter, xout) # TODO: 3rd order only!
                    push!(tspan, tspan_fixed[iter_fixed])
                    iter_fixed += 1
                    iter += 1
                end
                # but also output every step taken
                index_or_push!(xs, iter, copy(xtrial))
                push!(tspan, t+dt)
                iter += 1
            end
            # k[1] = f1 # load k[1]==f0 for next step

            # Break if this was the last step:
            laststep && break

            # Swap bindings of x and xtrial, avoids one copy
            x, xtrial = xtrial, x

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
    return tspan, xs

end
function rosw_step!{N,S}(xtrial, g!, gprime!, x, dt, dof, btab::TableauRosW_T{N,S},
                         k, jac_store, ks, u, udot, bt_ind=1)
    # This takes one step for a ode/dae system defined by
    # g(x,xdot)=0
    # gprime(x, xdot, α) = dg/dx + α dg/dxdot

    ## Calculate k
    # first step:
    s = 1
    h!(k[s], ks, u, udot, g!, x, btab, dt, k, s, dof) # k[s] now holds residual
    # It's sufficient to only update the Jacobian once per time step:
    # (Note that this uses u & udot calculated in h! above.)
    jac = hprime!(jac_store, x, gprime!, u, udot, btab, dt)
    # middle steps:
    for s=1:S-1
        # first step of Newton iteration with guess ks==0
        #        k[s][:] = ks - jac\k[s] # in-place A_ldiv_B!(jac, k[s])
        # TODO: add option to do more Newton iterations
        A_ldiv_B!(jac, k[s])
        for d=1:dof
            k[s][d] = ks[d]-k[s][d]
        end
        h!(k[s+1], ks, u, udot, g!, x, btab, dt, k, s+1, dof)
        # TODO: add option to evaluate Jakobian again
    end
    # last step:
    #    k[S][:] = ks - jac\k[S]
    # TODO: add option to do more Newton iterations
    A_ldiv_B!(jac, k[S])
    for d=1:dof
        k[S][d] = ks[d]-k[S][d]
    end
    
    ## Completion:
    xtrial[:] = x
    for s=1:S
        for d=1:dof
            xtrial[d] += btab.b[bt_ind,s]*k[s][d]
        end
    end
    return nothing
end
function h!(res, ks, u, udot, g!, x, btab, dt, k, s, dof)
    # h(ks)=0 to be solved for ks.
    #
    # Calculates h(ks) and h'(ks) (if s==1).
    # Modifies its first 3 arguments.
    #
    # res -- result: can use k[s] for this
    # ks -- guess for k[s] (usually==0) AND result g(u,udot)
    # u, udot -- work arrays for first and second argument to g
    # g -- obj function
    # x -- current state
    # btab
    # k -- stage vars (jed's yi)
    # s -- current stage to calculate
    # dof

    # stage independent
    for d=1:dof
        u[d] = ks[d] + x[d]
        udot[d] = 1/(dt*btab.γii)*ks[d] # usually ==0 as ks==0
    end
    # stages
    for ss=1:s-1
        for d=1:dof
            u[d] += btab.a[s,ss]*k[ss][d]
            udot[d] += btab.c[s,ss]/dt*k[ss][d]
        end
    end
    g!(res, u, udot)
    return nothing
end
function hprime!(res, xi, gprime!, u, udot, btab, dt)
    # The Jacobian of h!  Note that u and udot need to be calculated
    # by running h! first.
    #
    # The res will hold part of the LU factorization.  However use the
    # returned LU-type.
    gprime!(res, u, udot, 1./(dt*btab.γii))
    return lufact!(res)
end

# Copied & adapted from DASSL:
# generate a function that computes approximate jacobian using forward
# finite differences
function numerical_jacobian(fn!, reltol, abstol)
    function numjac!(jac, y, dy, a)
        # TODO use coloring
        
        ep      = eps(1e6)   # this is the machine epsilon # TODO: better value?
        h       = 1/a           # h ~ 1/a
        wt      = reltol*abs(y).+abstol
        # delta for approximation of jacobian.  I removed the
        # sign(h_next*dy0) from the definition of delta because it was
        # causing trouble when dy0==0 (which happens for ord==1)
        edelta  = spdiagm(max(abs(y),abs(h*dy),wt)*sqrt(ep))

        tmp = similar(y)
        f0 = similar(y)
        fn!(f0, y, dy)
        for i=1:length(y)
            fn!(tmp, y+edelta[:,i], a*edelta[:,i]+dy)
            jac[:,i] = (tmp-f0)/edelta[i,i]
        end
        return nothing
    end
end
