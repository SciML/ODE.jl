# Explicit Runge-Kutta solvers
##############################
# (Hairer & Wanner 1992 p.134, p.165-169)

###########################################
# Tableaus for explicit Runge-Kutta methods
###########################################

immutable TableauRKExplicit{Name, S, T} <: Tableau{Name, S, T}
    order::(@compat(Tuple{Vararg{Int}})) # the order of the methods
    a::Matrix{T}
    # one or several row vectors.  First row is used for the step,
    # second for error calc.
    b::Matrix{T}
    c::Vector{T}
    function TableauRKExplicit(order,a,b,c)
        @assert isa(S,Integer)
        @assert isa(Name,Symbol)
        @assert c[1]==0
        @assert istril(a)
        @assert S==length(c)==size(a,1)==size(a,2)==size(b,2)
        @assert size(b,1)==length(order)
        @assert norm(sum(a,2)-c'',Inf)<1e-10 # consistency.  
        new(order,a,b,c)
    end
end
function TableauRKExplicit{T}(name::Symbol, order::(@compat(Tuple{Vararg{Int}})),
                   a::Matrix{T}, b::Matrix{T}, c::Vector{T})
    TableauRKExplicit{name,length(c),T}(order, a, b, c)
end
function TableauRKExplicit(name::Symbol, order::(@compat(Tuple{Vararg{Int}})), T::Type,
                   a::Matrix, b::Matrix, c::Vector)
    TableauRKExplicit{name,length(c),T}(order, convert(Matrix{T},a),
                                        convert(Matrix{T},b), convert(Vector{T},c) )
end
conv_field{T,N}(D,a::Array{T,N}) = convert(Array{D,N}, a)
function Base.convert{Tnew<:Real,Name,S,T}(::Type{Tnew}, tab::TableauRKExplicit{Name,S,T})
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
    TableauRKExplicit{Name,S,Tnew}(newflds...) # TODO: could this be done more generically in a type-stable way?
end


isexplicit(b::TableauRKExplicit) = istril(b.a) # Test whether it's an explicit method
isadaptive(b::TableauRKExplicit) = size(b.b, 1)==2


# First same as last.  Means ks[:,end]=ks_nextstep[:,1], c.f. H&W p.167
isFSAL(btab::TableauRKExplicit) = btab.a[end,:]==btab.b[1,:] && btab.c[end]==1 # the latter is not needed really

## Tableaus for explicit RK methods
# Fixed step:
const bt_feuler = TableauRKExplicit(:feuler,(1,), Rational{Int64},
                                    zeros(Int,1,1),
                                    [1]',
                                    [0]
                                    )
const bt_midpoint = TableauRKExplicit(:midpoint,(2,), Rational{Int64},
                                      [0  0
                                       1//2  0],
                                      [0, 1]',
                                      [0, 1//2]
                                      )
const bt_heun = TableauRKExplicit(:heun,(2,), Rational{Int64},
                                  [0  0
                                   1  0],
                                  [1//2, 1//2]',
                                  [0, 1])

const bt_rk4 = TableauRKExplicit(:rk4,(4,),Rational{Int64},
                                 [0    0    0 0
                                  1//2 0    0 0
                                  0    1//2 0 0
                                  0    0    1 0],
                                 [1//6, 1//3, 1//3, 1//6]',
                                 [0, 1//2, 1//2, 1])

# Adaptive step:
# Heun Euler https://en.wikipedia.org/wiki/Runge–Kutta_methods
const bt_rk21 = TableauRKExplicit(:heun_euler,(2,1), Rational{Int64},
                                  [0     0
                                   1     0],
                                  [1//2  1//2
                                   1     0],
                                  [0,    1])

# Bogacki–Shampine coefficients
const bt_rk23 = TableauRKExplicit(:bogacki_shampine,(2,3), Rational{Int64},
                                  [0           0      0      0
                                   1/2         0      0      0
                                   0         3/4      0      0
                                   2/9       1/3     4/9     0],
                                  [7/24 1/4 1/3 1/8
                                   2/9 1/3 4/9 0],
                                  [0, 1//2, 3//4, 1] 
                         )

# Fehlberg https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method
const bt_rk45 = TableauRKExplicit(:fehlberg,(4,5),Rational{Int64},
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
const bt_dopri5 = TableauRKExplicit(:dopri, (5,4), Rational{Int64},
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
const bt_feh78 = TableauRKExplicit(:feh78, (7,8), Rational{Int64},
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


################################
# Fixed step Runge-Kutta methods
################################

# TODO: iterator method
ode1(fn, y0, tspan) = oderk_fixed(fn, y0, tspan, bt_feuler)
ode2_midpoint(fn, y0, tspan) = oderk_fixed(fn, y0, tspan, bt_midpoint)
ode2_heun(fn, y0, tspan) = oderk_fixed(fn, y0, tspan, bt_heun)
ode4(fn, y0, tspan) = oderk_fixed(fn, y0, tspan, bt_rk4)

function oderk_fixed(fn, y0, tspan, btab::TableauRKExplicit)
    # Non-arrays y0 treat as scalar
    fn_(t, y) = [fn(t, y[1])]
    t,y = oderk_fixed(fn_, [y0], tspan, btab)
    return t, vcat_nosplat(y)
end
function oderk_fixed{N,S}(fn, y0::AbstractVector, tspan,
                          btab_::TableauRKExplicit{N,S})
    # TODO: instead of AbstractVector use a Holy-trait

    # Needed interface:
    # On components: 
    # On y0 container: length, deepcopy, similar, setindex!
    # On time container: getindex, convert. length

    Et, Eyf, Ty, btab = make_consistent_types(fn, y0, tspan, btab_)
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
        ys[i+1][:] = ys[i]
        for s=1:S
            calc_next_k!(ks, ytmp, ys[i], s, fn, tspan[i], dt, dof, btab)
            for d=1:dof
                ys[i+1][d] += dt * btab.b[s]*ks[s][d]
            end
        end
    end
    return tspan, ys
end

##############################
# Adaptive Runge-Kutta methods
##############################

ode21(fn, y0, tspan; kwargs...) = oderk_adapt(fn, y0, tspan, bt_rk21; kwargs...)
ode23(fn, y0, tspan; kwargs...) = oderk_adapt(fn, y0, tspan, bt_rk23; kwargs...)
ode45_fe(fn, y0, tspan; kwargs...) = oderk_adapt(fn, y0, tspan, bt_rk45; kwargs...)
ode45_dp(fn, y0, tspan; kwargs...) = oderk_adapt(fn, y0, tspan, bt_dopri5; kwargs...)
# Use Dormand-Prince version of ode45 by default
const ode45 = ode45_dp
ode78(fn, y0, tspan; kwargs...) = oderk_adapt(fn, y0, tspan, bt_feh78; kwargs...)

function oderk_adapt(fn, y0, tspan, btab::TableauRKExplicit; kwords...)
    # For y0 which don't support indexing.
    fn_ = (t, y) -> [fn(t, y[1])]
    t,y = oderk_adapt(fn_, [y0], tspan, btab; kwords...)
    return t, vcat_nosplat(y)
end
function oderk_adapt{N,S}(fn, y0::AbstractVector, tspan, btab_::TableauRKExplicit{N,S};
                          reltol = 1.0e-5, abstol = 1.0e-8,
                          norm=Base.norm,
                          minstep=abs(tspan[end] - tspan[1])/1e18,
                          maxstep=abs(tspan[end] - tspan[1])/2.5,
                          initstep=0.,
                          points=:all
                          )
    # Needed interface:
    # On components:
    #  - note that the type of the components might change!
    # On y0 container: length, similar, setindex!
    # On time container: getindex, convert, length
    
    # For y0 which support indexing.  Currently y0<:AbstractVector but
    # that could be relaxed with a Holy-trait.
    !isadaptive(btab_) && error("Can only use this solver with an adaptive RK Butcher table")

    Et, Eyf, Ty, btab = make_consistent_types(fn, y0, tspan, btab_)
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
    dt, tdir, ks[1] = hinit(fn, y, tstart, tend, order, reltol, abstol) # sets ks[1]=f0
    if initstep!=0
        dt = sign(initstep)==tdir ? initstep : error("initstep has wrong sign.")
    end
    # Diagnostics
    dts = Et[]
    errs = Float64[]
    steps = [0,0]  # [accepted, rejected]

    ## Integration loop
    islaststep = abs(t+dt-tend)<=eps(tend) ? true : false
    timeout = 0 # for step-control
    iter = 2 # the index into tspan and ys
    while true
        # do one step (assumes ks[1]==f0)
        rk_embedded_step!(ytrial, yerr, ks, ytmp, y, fn, t, dt, dof, btab)
        # Check error and find a new step size:
        err, newdt, timeout = stepsize_hw92!(dt, tdir, y, ytrial, yerr, order, timeout,
                                            dof, abstol, reltol, maxstep, norm)

        if err<=1.0 # accept step
            # diagnostics
            steps[1] +=1
            push!(dts, dt)
            push!(errs, err)

            # Output:
            f0 = ks[1]
            f1 = isFSAL(btab) ? ks[S] : fn(t+dt, ytrial)
            if points==:specified
                # interpolate onto given output points
                while iter-1<nsteps_fixed && (tdir*tspan[iter]<tdir*(t+dt) || islaststep) # output at all new times which are < t+dt
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
                index_or_push!(ys, iter, deepcopy(ytrial))
                push!(tspan, t+dt)
                iter += 1
            end
            ks[1] = f1 # load ks[1]==f0 for next step

            # Break if this was the last step:
            islaststep && break

            # Swap bindings of y and ytrial, avoids one copy
            y, ytrial = ytrial, y

            # Update t to the time at the end of current step:
            t += dt
            dt = newdt

            # Hit end point exactly if next step within 1% of end:
            if tdir*(t+dt*1.01) >= tdir*tend
                dt = tend-t
                islaststep = true # next step is the last, if it succeeds
            end
        elseif abs(newdt)<minstep  # minimum step size reached, break
            println("Warning: dt < minstep.  Stopping.")
            break
        else # redo step with smaller dt
            islaststep = false
            steps[2] +=1
            dt = newdt
            timeout = timeout_const
        end
    end
    return tspan, ys
end

function rk_embedded_step!{N,S}(ytrial, yerr, ks, ytmp, y, fn, t, dt, dof, btab::TableauRKExplicit{N,S})
    # Does one embedded R-K step updating ytrial, yerr and ks.
    #
    # Assumes that ks[:,1] is already calculated!
    #
    # Modifies ytrial, yerr, ks, and ytmp

    # Needed interface:
    # On components: arithmetic, zero
    # On y0 container: fill!, setindex!, getindex

    fill!(ytrial, zero(eltype(ytrial)) )
    fill!(yerr, zero(eltype(ytrial)) )
    for d=1:dof
        ytrial[d] += btab.b[1,1]*ks[1][d]
        yerr[d]   += btab.b[2,1]*ks[1][d]
    end
    for s=2:S
        calc_next_k!(ks, ytmp, y, s, fn, t, dt, dof, btab)
        for d=1:dof
            ytrial[d] += btab.b[1,s]*ks[s][d]
            yerr[d]   += btab.b[2,s]*ks[s][d]
        end
    end
    for d=1:dof
        yerr[d]   = dt * (ytrial[d]-yerr[d])
        ytrial[d] = y[d] + dt * ytrial[d]
    end
end

function stepsize_hw92!(dt, tdir, x0, xtrial, xerr, order,
                       timeout, dof, abstol, reltol, maxstep, norm)
    # Estimates the error and a new step size following Hairer &
    # Wanner 1992, p167 (with some modifications)
    #
    # If timeout>0 no step size increase is allowed, timeout is
    # decremented in here.
    #
    # Returns the error, newdt and the number of timeout-steps
    #
    # TODO:
    # - allow component-wise reltol and abstol?
    # - allow other norms

    # Needed interface:
    # On components: isoutofdomain
    # On y0 container: norm, get/setindex

    timout_after_nan = 5
    fac = [0.8, 0.9, 0.25^(1/(order+1)), 0.38^(1/(order+1))][1]
    facmax = 5.0 # maximal step size increase. 1.5-5
    facmin = 1./facmax  # maximal step size decrease. ?

    # in-place calculate xerr./tol
    for d=1:dof
        # if outside of domain (usually NaN) then make step size smaller by maximum
        isoutofdomain(xtrial[d]) && return 10., dt*facmin, timout_after_nan
        xerr[d] = xerr[d]/(abstol + max(norm(x0[d]), norm(xtrial[d]))*reltol) # Eq 4.10
    end
    err = norm(xerr, 2) # Eq. 4.11
    newdt = min(maxstep, tdir*dt*max(facmin, fac*(1/err)^(1/(order+1)))) # Eq 4.13 modified
    if timeout>0
        newdt = min(newdt, dt)
        timeout -= 1
    end
    return err, tdir*newdt, timeout
end

function calc_next_k!{Ty}(ks::Vector, ytmp::Ty, y, s, fn, t, dt, dof, btab)
    # Calculates the next ks and puts it into ks[s]
    # - ks and ytmp are modified inside this function.

    # Needed interface:
    # On components: +, *
    # On y0 container: setindex!, getindex, fn

    ytmp[:] = y
    for ss=1:s-1, d=1:dof
        ytmp[d] += dt * ks[ss][d] * btab.a[s,ss]
    end
    ks[s] = fn(t + btab.c[s]*dt, ytmp)::Ty
    nothing
end

# Helper functions:
function allocate!{T}(vec::Vector{T}, y0, dof)
    # Allocates all vectors inside a Vector{Vector} using the same
    # kind of container as y0 has and element type eltype(eltype(vec)).
    for s=1:length(vec)
        vec[s] = similar(y0, eltype(T), dof)
    end
end
function index_or_push!(vec, i, val)
    # Fills in the vector until there is no space, then uses push!
    # instead.
    if length(vec)>=i
        vec[i] = val
    else
        push!(vec, val)
    end
end
vcat_nosplat(y) =  eltype(y[1])[el[1] for el in y] # Does vcat(y...) without the splatting

function hermite_interp!(y, tquery,t,dt,y0,y1,f0,f1)
    # For dense output see Hairer & Wanner p.190 using Hermite
    # interpolation. Updates y in-place.
    #
    # f_0 = f(x_0 , y_0) , f_1 = f(x_0 + h, y_1 )
    # this is O(3). TODO for higher order.

    theta = (tquery-t)/dt
    for i=1:length(y0)
        y[i] = ((1-theta)*y0[i] + theta*y1[i] + theta*(theta-1) *
                ((1-2*theta)*(y1[i]-y0[i]) + (theta-1)*dt*f0[i] + theta*dt*f1[i]) )
    end
    nothing
end
function hermite_interp(tquery,t,dt,y0,y1,f0,f1)
    # Returns the y instead of in-place
    y = similar(y0)
    hermite_interp!(y,tquery,t,dt,y0,y1,f0,f1)
    return y
end
