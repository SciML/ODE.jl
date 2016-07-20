abstract Options{T}

"""

Generic options for adaptive ODE solvers.  This type has a key-word
constructor which will fill the structure with default values.

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
immutable AdaptiveOptions{T,N<:Function,O<:Function} <: Options{T}
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

@compat function (::Type{AdaptiveOptions{T}}){T,N,O}(;
                                                     tspan    = T[Inf],
                                                     tstop    = tspan[end],
                                                     reltol   = eps(T)^T(1//3)/10,
                                                     abstol   = eps(T)^T(1//2)/10,
                                                     minstep  = 10*eps(T),
                                                     maxstep  = 1/minstep,
                                                     initstep = minstep,
                                                     norm::N  = Base.norm,
                                                     maxiters = T(Inf),
                                                     isoutofdomain::O = Base.isnan,
                                                     kargs...)

    AdaptiveOptions{T,N,O}(tstop,reltol,abstol,minstep,maxstep,initstep,norm,maxiters,isoutofdomain)
end

"""

Generic options for fixed step ODE solvers.  This type has a key-word
constructor which will fill the structure with default values.

General:

- initstep ::T  initial step
- tstop    ::T  end integration time

"""
immutable FixedOptions{T} <: Options{T}
    tstop::T
    initstep::T
end

@compat function (::Type{FixedOptions{T}}){T}(;
                                              tspan    = T[Inf],
                                              tstop    = tspan[end],
                                              initstep = 10*eps(T),
                                              kargs...)

    FixedOptions{T}(tstop,initstep)
end

function show{T}(io::IO, opts :: Options{T})
    for name in fieldnames(opts)
        @printf("%-20s = %s\n",name,getfield(opts,name))
    end
end
