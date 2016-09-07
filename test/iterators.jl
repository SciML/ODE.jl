
const integrators = [ODE.RKIntegratorFixed{:feuler},
                  ODE.RKIntegratorFixed{:midpoint},
                  ODE.RKIntegratorFixed{:heun},
                  ODE.RKIntegratorFixed{:rk4},
                  ODE.RKIntegratorAdaptive{:rk21},
                  ODE.RKIntegratorAdaptive{:rk23},
                  ODE.RKIntegratorAdaptive{:rk45},
                  ODE.RKIntegratorAdaptive{:dopri5},
                  ODE.RKIntegratorAdaptive{:feh78},
                  ODE.ModifiedRosenbrockIntegrator
                  ]

@testset "Iterator interfaces" begin
    for integ in integrators
        ODETests.test_integrator(integ)
    end
end
