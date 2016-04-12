const steppers =
    [
     ( :ode23s,        :ModifiedRosenbrockStepper, []),
     ( :ode1,          :TableauStepperFixed,       [bt_feuler]),
     ( :ode2_midpoint, :TableauStepperFixed,       [bt_midpoint]),
     ( :ode2_heun,     :TableauStepperFixed,       [bt_heun]),
     ( :ode4,          :TableauStepperFixed,       [bt_rk4]),
     ( :ode21,         :TableauStepperAdaptive,    [bt_rk21]),
     ( :ode23,         :TableauStepperAdaptive,    [bt_rk23]),
     ( :ode45_fe,      :TableauStepperAdaptive,    [bt_rk45]),
     ( :ode45_dp,      :TableauStepperAdaptive,    [bt_dopri5]),
     ( :ode78,         :TableauStepperAdaptive,    [bt_feh78])
]



# TODO: there is a lot of useless conversions going on here

for (name,stepper,params) in steppers
    @eval begin
        function ($name){T<:Number}(F, y0 :: Vector, t0 :: T;
                                    jacobian = (t,y)->fdjacobian(F, t, y),
                                    stopevent = (t,y)->false,
                                    tstop = Inf,
                                    tspan = [tstop],
                                    kargs...)

            step = ($stepper){T}($params...)

            if all(tspan .>= t0)
                # forward time integration
                ode  = explicit_ineff(t0,y0,F,jac=jacobian)
                opts = Options{T}(;
                                  tstop = tstop,
                                  tspan = tspan,
                                  stopevent = stopevent,
                                  kargs...)
                solution = collect(dense(solve(ode,step,opts)))
                n  = length(solution)
            elseif all(tspan .<= t0)
                # reverse time integration
                F_reverse(t,y) = -F(2*t0-t,y)
                # TODO: is that how the jacobian changes?
                jac_reverse(t,y) = -jacobian(2*t0-t,y)
                ode  = explicit_ineff(t0,y0,F_reverse,jac=jac_reverse)
                opts = Options{T}(;
                                  tstop = 2*t0-tstop,
                                  tspan = 2*t0.-tspan,
                                  stopevent = (t,y)->stopevent(2*t0-t,y),
                                  kargs...)
            else
                # tspan stretches to the left and to the right of t0
                return ([t0],[y0])
            end

            solution = collect(dense(solve(ode,step,opts)))
            n = length(solution)

            # return solution

            # convert a list of pairs to a pair of arrays
            # TODO: leave it out as a list of pairs?
            tn = Array(T,n)
            yn = Array(typeof(y0),n)
            if all(tspan .>= t0)
                tn[:] = [x[1] for x in solution]
            else
                tn[:] = [2*t0-x[1] for x in solution]
            end
            yn[:] = [x[2] for x in solution]

            return (tn,yn)
        end

        ($name){T<:Number}(F, y0 :: Vector, t0 :: Vector{T}; kargs...) =
            ($name)(F,y0,t0[1];
                    tstop  = t0[end],
                    tspan  = t0,
                    points = :specified,
                    kargs...)

        function ($name)(F, y0, t0; kargs...)
            tn, yn = ($name)((t,y)->[F(t,y[1])], [y0], t0; kargs...)
            yn2  = Array(typeof(y0),length(yn))
            yn2[:] = map(first,yn)
            return (tn,yn2)
        end
    end
end


const ode45 = ode45_dp
