function solve(
    prob::AbstractODEProblem{uType,tType,isinplace},
    alg::AlgType,
    timeseries=[], ts=[], ks=[];

    verbose=true,
    save_timeseries=nothing,
    saveat=tType[], reltol = 1e-5, abstol = 1e-8,
    save_everystep=isempty(saveat), 
    dense = save_everystep && !(typeof(alg) <: FunctionMap) && isempty(saveat),
    save_start = save_everystep || isempty(saveat) || typeof(saveat) <: Number ? true : prob.tspan[1] in saveat,
    callback=nothing,
    dtmin = abs(prob.tspan[2]-prob.tspan[1])/1e-9,
    dtmax = abs(prob.tspan[2]-prob.tspan[1])/2.5,
    timeseries_errors=true, dense_errors=false,
    dt = 0.0, norm = Base.vecnorm,
    kwargs...) where {uType,tType,isinplace,AlgType<:ODEjlAlgorithm}

    if verbose
        warned = !isempty(kwargs) && check_keywords(alg, kwargs, warnlist)
        if !(typeof(prob.f) <: AbstractParameterizedFunction) && typeof(alg) <: ode23s
            if has_tgrad(prob.f)
                warn("Explicit t-gradient given to this stiff solver is ignored.")
                warned = true
            end
            if has_jac(prob.f)
                warn("Explicit Jacobian given to this stiff solver is ignored.")
                warned = true
            end
        end
        warned && warn_compat()
    end

    if save_timeseries != nothing
        verbose && warn("save_timeseries is deprecated. Use save_everystep instead")
        save_everystep = save_timeseries
    end

    if prob.mass_matrix != I
        error("This solver is not able to use mass matrices.")
    end

    if callback != nothing || prob.callback != nothing
        error("ODE is not compatible with callbacks.")
    end

    tspan = prob.tspan

    u0 = prob.u0

    if typeof(saveat) <: Number
      if (tspan[1]:saveat:tspan[end])[end] == tspan[end]
        saveat_vec = convert(Vector{tType},collect(tType,tspan[1]+saveat:saveat:tspan[end]))
      else
        saveat_vec = convert(Vector{tType},collect(tType,tspan[1]+saveat:saveat:(tspan[end]-saveat)))
      end
    else
        saveat_vec = convert(Vector{tType}, collect(saveat))
    end

    if !isempty(saveat_vec) && saveat_vec[end] == tspan[2]
        pop!(saveat_vec)
    end

    if !isempty(saveat_vec) && saveat_vec[1] == tspan[1]
        Ts = unique([saveat_vec;tspan[2]])
    else
        Ts = unique([tspan[1];saveat_vec;tspan[2]])
    end

    if save_everystep
        points = :all
    else
        points = :specified
    end

    sizeu = size(prob.u0)
    p = prob.p

    if isinplace
        f = (t,u) -> (du = zeros(u); prob.f(du,u,p,t); vec(du))
    elseif uType <: AbstractArray
        f = (t,u) -> vec(prob.f(reshape(u,sizeu),p,t))
    else
        f = (t,u) -> prob.f(u,p,t)
    end

    u0 = uType <: AbstractArray ? vec(prob.u0) : prob.u0

    # Calling the solver, i.e. if the algorithm is ode45,
    # then AlgType(...) is ode45(...)
    ts, timeseries_tmp = AlgType(f,u0,Ts;
                                 norm = norm,
                                 abstol=abstol,
                                 reltol=reltol,
                                 maxstep=dtmax,
                                 minstep=dtmin,
                                 initstep=dt,
                                 points=points)

    if save_start
        start_idx = 1 # The index to start making the timeseries from
    else
        start_idx = 2
        ts = ts[2:end]
    end

    # Reshape the result if needed
    if uType <: AbstractArray
        timeseries = Vector{uType}(0)
        for i=start_idx:length(timeseries_tmp)
            push!(timeseries,reshape(timeseries_tmp[i],sizeu))
        end
    else
        timeseries = timeseries_tmp
    end

    build_solution(prob,alg,ts,timeseries,
                   timeseries_errors = timeseries_errors,
                   dense_errors = dense_errors,
                   retcode = :Succss)
end # function solve
