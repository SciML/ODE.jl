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

@compat function (::Type{StepperOptions{T}}){T,N,O}(;
                                                    tstop    = T(Inf),
                                                    reltol   = eps(T)^T(1//3)/10,
                                                    abstol   = eps(T)^T(1//2)/10,
                                                    minstep  = 10*eps(T),
                                                    maxstep  = 1/minstep,
                                                    initstep = minstep,
                                                    norm::N  = Base.norm,
                                                    steps    = repeated(initstep),
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
