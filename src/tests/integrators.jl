using Base.Test

"""
    test_integrator(integrator, [test_case])

Test a single `integrator` with either a single `test_case` or all
test cases.  Test cases are defined in `ODETests.test_cases` in
`src/tests/test_cases.jl`.

"""

function test_integrator(integrator,test)

    ivp, sol, name, opts = test[:ivp], test[:sol], test[:name], test[:options]

    T,Y = eltype(ivp).parameters

    tol = 1//500
    @testset "$name" begin

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
        iterator       = ODE.iterate(ivp; solver=integrator, opts...)
        iterator_dense = ODE.iterate(ivp; solver=ODE.DenseOutput{integrator}, opts...)

        tdir = sign(opts[:tstop]-ivp.t0)

        for iter in (iterator,iterator_dense)
            for (t,y,dy) in iter
                # TODO: is there a better way of doing this?  We need
                # a single @test statement in case of a failure for
                # @testset to work properly.
                if maxabs(y-sol(t)) > tol
                    @test maxabs(y-sol(t)) <= tol
                    break
                end
                if tdir*t > tdir*opts[:tstop]
                    @test tdir*t <= tdir*opts[:tstop]
                    break
                end
            end

            # generator comprehension
            @test all(collect((maxabs(y-sol(t))<=tol for (t,y) in iter)))
        end

        tout = opts[:tout]
        @test collect((t for (t,y) in iterator_dense))==tout

        # Solution type
        solution       = ODE.solve(ivp; solver=integrator, opts...)
        solution_dense = ODE.solve(ivp; solver=ODE.DenseOutput{integrator}, opts...)

        for s in (solution,solution_dense)
            @test all(map((t,y)->maxabs(y-sol(t))<=tol,solution.t,solution.y))
        end
    end
end

function test_integrator(integrator)
    @testset "$integrator" begin
        for case in values(test_cases)
            ODETests.test_integrator(integrator,case)
        end
    end
end

# this is the functionality not yet included in test_integrator
function unused()
# 3) test the backend API
    if properties != nothing
        order, name, isadaptive = properties
        @test ODE.order(integ) == order
        @test ODE.name(integ) == name
        @test ODE.isadaptive(integ) == isadaptive
    end
end
