using Base.Test

function test_integrator(integrator,test)

    ivp, sol, name, opts = test[:ivp], test[:sol], test[:name], test[:options]

    println("Integrator $integrator)")
    print("   Test case $name   ")

    T,Y = eltype(ivp).parameters

    tol = 1//500

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
    iterator = ODE.solve(ivp,integrator; opts...)
    dense    = ODE.solve(ivp,ODE.DenseOutput{integrator}; opts...)
    niters=0
    for (t,y,dy) in iterator
        niters+=1
        @test maxabs(y-sol(t)) < niters*tol
        # TODO: replace with
        # @test_approx_eq_eps y sol(t) tol
    end

    # with dense output
    for (t,y,dy) in dense
        @test maxabs(y-sol(t)) < tol
    end

    # generator comprehension
    @test all(collect((maxabs(y-sol(t))<=tol for (t,y) in iterator)))
    @test all(collect((maxabs(y-sol(t))<=tol for (t,y) in dense)))

    tout = opts[:tout]
    @test collect((t for (t,y) in dense))==tout

    println("OK!")
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
