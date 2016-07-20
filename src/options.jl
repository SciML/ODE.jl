abstract Options{T}

#m3: these are default options for adaptive stepper
"""

Options for ODE solvers.  This type has a key-word constructor which
will fill the structure with default values.

General:

- initstep ::T  initial step
- tstop    ::T  end integration time
- reltol   ::T  relative tolerance (m3: could this be a vector?)
- abstol   ::T  absolute tolerance (m3: could this be a vector?)
- minstep  ::T  minimal allowed step
- maxstep  ::T  maximal allowed step
- norm           function to calculate the norm in step control
- maxiters ::T  maximum number of steps
- isoutofdomain::Function checks if the solution is outside of the allowed domain

"""
immutable StepperOptions{T<:Number,N<:Function,O<:Function} <: Options{T}
    tstop::T
    reltol::T
    abstol::T
    minstep::T
    maxstep::T
    initstep::T
    norm::N
    maxiters::T
    isoutofdomain::O
end

@compat function (::Type{StepperOptions{T}}){T,N,O}(ode::ExplicitODE,
                                                    order::Int;
                                                    tspan::Vector = T[Inf],
                                                    tstop    = tspan[end],
                                                    reltol   = eps(T)^T(1//3)/10,
                                                    abstol   = eps(T)^T(1//2)/10,
                                                    minstep  = 10*eps(T),
                                                    maxstep  = 1/minstep,
                                                    initstep = dtinit(ode,order,reltol,abstol,tstop),
                                                    norm::N  = Base.norm,
                                                    maxiters = T(Inf),
                                                    isoutofdomain::O = Base.isnan,
                                                    kargs...)

    StepperOptions{T,N,O}(tstop,reltol,abstol,minstep,maxstep,initstep,norm,maxiters,isoutofdomain)
end

function show{T}(io::IO, opts :: Options{T})
    for name in fieldnames(opts)
        @printf("%-20s = %s\n",name,getfield(opts,name))
    end
end

#pwl: should this be dropped?
function dtinit{T}(ode::ExplicitODE{T},order::Int,reltol::T,abstol::T,tstop::T)
    t0 = ode.t0
    y0 = ode.y0

    f0 = similar(y0)
    tau = max(reltol*norm(y0, Inf), abstol)
    d0 = norm(y0, Inf)/tau
    ode.F!(t0, y0, f0)
    d1 = norm(f0, Inf)/tau
    if min(d0,d1) < eps(T)^(1/3)
        dt0 = eps(T)^(1/3)/10
    else
        dt0 = (d0/d1)/100
    end
    # perform Euler step
    y1 = y0+dt0*f0
    f1 = similar(f0)
    ode.F!(t0 + dt0, y1, f1)
    # estimate second derivative
    d2 = norm(f1 - f0, Inf)/(tau*dt0)
    if max(d1, d2) <= 10*eps(T)
        dt1 = max(eps(T)^(1/3)/10, dt0/10^3)
    else
        pow = -(2 + log10(max(d1, d2)))/(order+1)
        dt1 = 10^pow
    end
    return T(min(100*dt0, dt1, abs(tstop-t0)))
end
