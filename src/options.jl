abstract Options{T}

"""

Generic options for adaptive ODE solvers.  This type has a key-word
constructor which will fill the structure with default values.

General:

- initstep ::T  initial step size (always positive)
- tstop    ::T  end integration time
- reltol   ::T  relative tolerance (m3: could this be a vector?)
- abstol   ::T  absolute tolerance (m3: could this be a vector?)
- minstep  ::T  minimal allowed step size (always positive)
- maxstep  ::T  maximal allowed step size (always positive)
- norm          function to calculate the norm in step control
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
                                                     tout     = [T(1//0)],
                                                     tstop    = tout[end],
                                                     reltol   = eps(T)^T(1//3)/T(10),
                                                     abstol   = eps(T)^T(1//2)/T(10),
                                                     minstep  = T(10)*eps(T),
                                                     maxstep  = 1/minstep,
                                                     initstep = eps(T)^T(1//3),
                                                     norm::N  = y->maxabs(y),
                                                     maxiters = T(1//0),
                                                     isoutofdomain::O = Base.isnan,
                                                     kargs...)
    @assert minstep>=T(0) && maxstep>=T(0) && initstep>=T(0) # TODO: move to inner constructor
    AdaptiveOptions{T,N,O}(tstop,reltol,abstol,minstep,maxstep,initstep,norm,maxiters,isoutofdomain)
end

"""

Generic options for fixed step ODE solvers.  This type has a key-word
constructor which will fill the structure with default values.

General:

- initstep ::T  initial step (always positive)
- tstop    ::T  end integration time

"""
immutable FixedOptions{T} <: Options{T}
    tstop::T
    initstep::T
end

@compat function (::Type{FixedOptions{T}}){T}(;
                                              tout     = [T(1//0)],
                                              tstop    = tout[end],
                                              initstep = T(1//100),
                                              kargs...)
    @assert initstep>=0
    FixedOptions{T}(tstop,initstep)
end

function show{T}(io::IO, opts :: Options{T})
    for name in fieldnames(opts)
        @printf("%-20s = %s\n",name,getfield(opts,name))
    end
end
