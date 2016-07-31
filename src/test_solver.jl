import Base.Test: @test, @test_approx_eq_eps

const testset = [
                 Dict(:ivp      => ExplicitODE(0.0,
                                               [0.0],
                                               (t,y,dy)->dy[1]=6.0,
                                               J! =(t,y,dy)->dy[1]=0.0),
                      :sol      => t->[6t],
                      :isscalar => true,
                      :name     => "y'=6",
                      :options  => Dict(:initstep => 0.1,
                                        :tstop => 1.0)
                      )
                 ]


function test_integrator{I<:AbstractIntegrator}(integrator::Type{I};
                                                params = nothing,
                                                eqntypes = [ExplicitODE])

    const tol = 0.02

    # filter the tests
    tests = filter(x->in(typeof(x[:ivp]),eqntypes),testset)

    for test in testset
        println("Testing problem $(test[:name])")

        ivp, opts, sol = test[:ivp], test[:options], test[:sol]

        # 1) test the constructor
        @test integrator <: AbstractIntegrator
        integ=integrator(ivp;opts...)
        @test typeof(integ)<:integrator

        # 2) test if the minimal backend is implemented
        state=init(ivp,integ)
        @test typeof(state)<:AbstractState
        @test onestep!(ivp,integ,state) == cont
        # TODO: we should implement outputing the initial data as the
        # first step
        # @test output(state) == (ode.t0,ode.y0,ode.dy0)
        @test typeof(output(state)) == eltype(ivp)

        # 3) test the info methods
        if params != nothing
            order, name, isadaptive = params
            @test ODE.order(integ) == order
            @test ODE.name(integ) == name
            @test ODE.isadaptive(integ) == isadaptive
        end

        # 4) test the iterator interface
        # pure integrator
        for (t,y,dy) in ODE.solve(ivp,integrator;opts...)
            @test_approx_eq_eps y sol(t) tol
        end
        # with dense output
        for (t,y) in ODE.solve(ivp,ODE.DenseOutput{integrator}; opts...)
            @test_approx_eq_eps y sol(t) tol
        end
    end
end
