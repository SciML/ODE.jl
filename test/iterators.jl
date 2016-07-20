testsets = [
            Dict(
                 :F!    => (t,y,dy)->dy[1]=6.0,
                 :y0    => [0.],
                 :tspan => [0:0.1:1;],
                 :jac   => (t,y,dy)->dy[1]=0.0,
                 :sol   => t->[6t],
                 :isscalar => true,
                 :name  => "y'=6t",
                 :initstep => 0.1),
            Dict(
                 :F!    => (t,y,dy)->dy[1]=2t,
                 :y0    => [0.],
                 :tspan => [0:0.001:1;],
                 :jac   => (t,y,dy)->dy[1]=0.0,
                 :sol   => t->[t^2],
                 :isscalar => true,
                 :name  => "y'=2t",
                 :initstep => 0.001),
            Dict(
                 :F!    => (t,y,dy)->dy[1]=y[1],
                 :y0    => [1.0],
                 :tspan => [0:0.001:1;],
                 :jac   => (t,y,dy)->dy[1]=1.0,
                 :sol   => t->[exp(t)],
                 :isscalar => true,
                 :name  => "y'=y",
                 :initstep => 0.001),
            Dict(
                 :F!    => (t,y,dy)->dy[1]=y[1],
                 :y0    => [1.0],
                 :tspan => [1:-0.001:0;],
                 :jac   => (t,y,dy)->dy[1]=1.0,
                 :sol   => t->[exp(t-1)],
                 :isscalar => true,
                 :name  => "y'=y backwards",
                 :initstep => 0.001),
            Dict(
                 :F!    => (t,y,dy)->(dy[1]=-y[2];dy[2]=y[1]),
                 :y0    => [1.0,2.0],
                 :tspan => [0:.1:1;],
                 :jac   => (t,y,dy)->copy!(dy,Float64[[0,1] [-1,0]]),
                 :sol   => t->[cos(t)-2*sin(t) 2*cos(t)+sin(t)],
                 :isscalar => false,
                 :name  => "pendulum",
                 :initstep => 0.001)
            ]


# Testing function ode
steppers = [ODE.RKStepperFixed{:feuler},
            ODE.RKStepperFixed{:midpoint},
            ODE.RKStepperFixed{:heun},
            ODE.RKStepperFixed{:rk4},
            ODE.RKStepperAdaptive{:rk21},
            ODE.RKStepperAdaptive{:rk23},
            ODE.RKStepperAdaptive{:rk45},
            ODE.RKStepperAdaptive{:dopri5},
            ODE.RKStepperAdaptive{:feh78},
            ODE.ModifiedRosenbrockStepper{}
]

function test_ode()
    tol = 0.002

    for rks in steppers
        println("Testing $rks")
        for ts in testsets
            println("Testing problem $(ts[:name])")

            tspan, h0, stepper = ts[:tspan], ts[:initstep], rks{eltype(ts[:tspan])}()

            y0, F!, jac!, sol = ts[:y0], ts[:F!], ts[:jac], ts[:sol]

            F(t,y) = (dy = similar(y); F!(t,y,dy); return dy)

            for points = [:specified, :all]
                if ts[:isscalar]
                    # test the ODE.odeXX scalar interface (if the equation is scalar)
                    Fscal = (t,y)->F(t,[y])[1]
                    y0scal = y0[1]
                    # with jacobian
                    tj,yj = ODE.ode(Fscal,y0scal,tspan,stepper,points=points,initstep = h0,J! = jac!)
                    @test_approx_eq_eps yj map(x->sol(x)[1],tj) tol
                    # without jacobian
                    t,y   = ODE.ode(Fscal,y0scal,tspan,stepper,points=points,initstep = h0)
                    @test_approx_eq_eps y  map(x->sol(x)[1],tj) tol

                    # results with and without jacobian should be exactly the same
                    @test_approx_eq yj y

                    if points == :specified
                        # test if we covered the whole timespan
                        @test length(tspan) == length(t) == length(tj)
                        @test_approx_eq tspan t
                        @test_approx_eq tspan tj
                    end
                end

                # ODE.odeXX vector interface
                # with jacobian
                tj,yj = ODE.ode(F,y0,tspan,stepper,points=points,initstep = h0,J! = jac!)
                @test_approx_eq_eps hcat(yj...) hcat(map(sol,tj)...) tol
                # without jacobian
                t,y   = ODE.ode(F,y0,tspan,stepper,points=points,initstep = h0)
                @test_approx_eq_eps hcat(y...)  hcat(map(sol,t)...) tol

                @test_approx_eq hcat(yj...) hcat(y...)

                if points == :specified
                    # test if we covered the whole timespan
                    @test length(tspan) == length(t) == length(tj)
                    @test_approx_eq tspan t
                    @test_approx_eq tspan tj
                end

                # test the iterator interface (they only support forward time integration)
                if issorted(tspan)
                    equation = ODE.ExplicitODE(tspan[1],y0,F!)
                    opts     = ODE.Options{eltype(tspan)}(tspan = tspan,initstep = h0,points = points)
                    solver   = ODE.solve(equation,stepper,opts)

                    for (t,y) in solver
                        @test_approx_eq_eps y sol(t) tol
                    end

                    for (t,y) in ODE.dense(solver)
                        @test_approx_eq_eps y sol(t) tol
                    end
                end
            end
        end
    end
end

test_ode()
