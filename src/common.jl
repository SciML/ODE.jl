function solve{uType,tType,isinplace,AlgType<:ODEJLAlgorithm,F}(prob::AbstractODEProblem{uType,tType,isinplace,F},
    alg::AlgType,timeseries=[],ts=[],ks=[];dense=true,save_timeseries=true,
    saveat=tType[],timeseries_errors=true,reltol = 1e-5, abstol = 1e-8,
    dtmin = abs(prob.tspan[2]-prob.tspan[1])/1e-9,
    dtmax = abs(prob.tspan[2]-prob.tspan[1])/2.5,
    dt = 0.,norm = Base.vecnorm,
    kwargs...)

    tspan = prob.tspan

    if tspan[end]-tspan[1]<tType(0)
        error("final time must be greater than starting time. Aborting.")
    end

    u0 = prob.u0

    Ts = sort(unique([tspan[1];saveat;tspan[2]]))

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

    if uType <: AbstractArray
        u0 = vec(prob.u0)
    else
        u0 = prob.u0
    end

    ts,timeseries_tmp = AlgType(f,u0,Ts;
                      norm = norm,
                      abstol=abstol,
                      reltol=reltol,
                      maxstep=dtmax,
                      minstep=dtmin,
                      initstep=dt,
                      points=points)
                      #=
    elseif typeof(alg) <: ode4
        ts,timeseries_tmp = ODE.ode4(f,u0,Ts)
    elseif typeof(alg) <: ode4ms
        ts,timeseries_tmp = ODE.ode4ms(f,u0,Ts)
    elseif typeof(alg) <: ode4s
        ts,timeseries_tmp = ODE.ode4s(f,u0,Ts)
    end
    =#

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
                 timeseries_errors = timeseries_errors)
end

export ODEJLAlgorithm, ode23Alg, ode23sAlg, ode45Alg, ode78Alg
