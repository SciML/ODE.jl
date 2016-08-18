using Base.Test

const case_vector = Dict(:ivp      => ExplicitODE(0.0,
                                                  [0.0],
                                                  (t,y,dy)->dy[:]=6.0,
                                                  J! =(t,y,dy)->dy[1]=0.0),
                         :sol      => t->[6t],
                         :name     => "y'=6 (vector)",
                         :options  => Dict(:initstep => 0.1,
                                           :tstop => 1.0),
                         )

function test_integrator(integrator,test)

    ivp, sol, name, opts = test[:ivp], test[:sol], test[:name], test[:options]

    T,Y = eltype(ivp).parameters

    tol = 2//10

    # 1) test the constructor
    @test integrator <: AbstractIntegrator
    integ=integrator(ivp;opts...)
    @test typeof(integ)<:integrator

    # 2) test if the minimal backend is implemented
    state=ODE.init(ivp,integ)
    @test typeof(state)<:AbstractState
    # output after initialization should give the initial data
    @test ODE.output(state) == (ivp.t0,ivp.y0,ivp.dy0)

    # we should be able to perform the first step
    @test ODE.onestep!(ivp,integ,state) == ODE.cont
    # after one step the output should be updated
    @test ODE.output(state) != (ivp.t0,ivp.y0,ivp.dy0)

    # 3) test the iterator interface
    # pure integrator
    for (t,y,dy) in ODE.solve(ivp,integrator;opts...)
        @test maxabs(y-sol(t)) < tol
    end
    # with dense output
    for (t,y) in ODE.solve(ivp,ODE.DenseOutput{integrator}; opts...)
        @test maxabs(y-sol(t)) < tol
        # TODO: make this work
        # @test_approx_eq_eps y sol(t) tol
    end
end

function aaaa()
# 3) test the backend API
    if properties != nothing
        order, name, isadaptive = properties
        @test ODE.order(integ) == order
        @test ODE.name(integ) == name
        @test ODE.isadaptive(integ) == isadaptive
    end
end
