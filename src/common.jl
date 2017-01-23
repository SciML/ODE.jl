function solve{uType,tType,isinplace,AlgType<:ODEjlAlgorithm,F}(prob::AbstractODEProblem{uType,tType,isinplace,F},
    alg::AlgType,timeseries=[],ts=[],ks=[];dense=true,save_timeseries=true,
    saveat=tType[],reltol = 1e-5, abstol = 1e-8,
    dtmin = abs(prob.tspan[2]-prob.tspan[1])/1e-9,
    dtmax = abs(prob.tspan[2]-prob.tspan[1])/2.5,
    timeseries_errors=true,dense_errors=false,
    dt = 0.,norm = Base.vecnorm,
    kwargs...)

    tspan = prob.tspan

    u0 = prob.u0

    Ts = unique([tspan[1];saveat;tspan[2]])

    if save_timeseries
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

    # Reshape the result if needed
    if uType <: AbstractArray
        timeseries = Vector{uType}(0)
        for i=1:length(timeseries_tmp)
            push!(timeseries,reshape(timeseries_tmp[i],sizeu))
        end
    else
        timeseries = timeseries_tmp
    end

    build_solution(prob,alg,ts,timeseries,
                 timeseries_errors = timeseries_errors,
                 dense_errors = dense_errors)
end
