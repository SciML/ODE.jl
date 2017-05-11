function solve{uType,tType,isinplace,AlgType<:ODEjlAlgorithm}(prob::AbstractODEProblem{uType,tType,isinplace},
    alg::AlgType,timeseries=[],ts=[],ks=[];dense=true,
    save_timeseries=nothing,
    saveat=tType[],reltol = 1e-5, abstol = 1e-8,
    save_everystep=isempty(saveat),
    save_start = true,callback=nothing,
    dtmin = abs(prob.tspan[2]-prob.tspan[1])/1e-9,
    dtmax = abs(prob.tspan[2]-prob.tspan[1])/2.5,
    timeseries_errors=true,dense_errors=false,
    dt = 0.,norm = Base.vecnorm,
    kwargs...)

    if save_timeseries != nothing
        warn("save_timeseries is deprecated. Use save_everystep instead")
        save_everystep = save_timeseries
    end
    
    if prob.mass_matrix != I
        error("This solver is not able to use mass matrices.")
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
      saveat_vec = convert(Vector{tType},saveat:saveat:(tspan[end]-saveat))
      # Exclude the endpoint because of floating point issues
    else
      saveat_vec =  convert(Vector{tType},collect(saveat))
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

    if isinplace
        f = (t,u) -> (du = zeros(u); prob.f(t,u,du); vec(du))
    elseif uType <: AbstractArray
        f = (t,u) -> vec(prob.f(t,reshape(u,sizeu)))
    else
        f = prob.f
    end

    u0 = uType <: AbstractArray ? vec(prob.u0) : prob.u0

    # Calling the solver, i.e. if the algorithm is ode45,
    # then AlgType(...) is ode45(...)
    ts,timeseries_tmp = AlgType(f,u0,Ts;
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
end
