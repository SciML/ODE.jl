# Test sets [F, y0, tspan, analytic]
testsets_scalar = Vector[
                         Any[(t,y)->6.0, 0., [0:.1:1;], (t,y)->6t],
                         Any[(t,y)->2t, 0., [0:.001:1;], (t,y)->t.^2],
                         Any[(t,y)->y, 1., [0:.001:1;], (t,y)->e.^t],
                         Any[(t,y)->y, 1., [1:-.001:0;], (t,y)->e.^(t-1)],
                         Any[(t,y)->[-y[2]; y[1]], [1., 2.], [0:.001:2*pi;],
                             (t,y)->[cos(t)-2*sin(t) 2*cos(t)+sin(t)] ]
]

testsets_vector = Vector[
                         Any[(t,y,dy)-> dy[1]=6.0, [0.], [0:.1:1;], (t,y)->6t],
                         Any[(t,y,dy)-> dy[1]=2t, [0.], [0:.001:1;], (t,y)->t.^2],
                         Any[(t,y,dy)-> dy[1]=y[1], [1.], [0:.001:1;], (t,y)->e.^t],
                         Any[(t,y,dy)-> dy[1]=y[1], [1.], [1:-.001:0;], (t,y)->e.^(t-1)],
                         Any[(t,y,dy)->(dy[1]=-y[2]; dy[2]=y[1]), [1., 2.], [0:.001:2*pi;],
                             (t,y)->[cos(t)-2*sin(t) 2*cos(t)+sin(t)] ]
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

# F,y0,tspan,ana = (1,1,1,1)
rks =1
ts =1
println("Testing `ode`")
function test_ode()
    for rks in steppers
        println("Testing $rks")
        for ts in testsets_scalar
            F,y0,tspan,ana = ts
            t,y = ODE.ode(F,y0,tspan,rks{eltype(tspan)}())
            y = hcat(y...).'
            @test maximum(abs(y-ana(t,y))) < tol
        end
    end
end

function test_iterator_out_place()
    # Testing the lower-level iteration API
    println("\nTesting iterators")
    for rks in steppers
        println("Testing $rks")
        for ts in testsets_scalar
            F,y0,tspan,ana = ts
            T = eltype(tspan)
            stepper = rks{T}()
            jac = ODE.forward_jacobian(F,y0)
            equation  = ODE.explicit_ineff(tspan[1],y0,F,jac)
            opts = ODE.Options{T}(
                                  tspan    = tspan,
                                  reltol   = eps(T)^T(1//3)/10,
                                  abstol   = eps(T)^T(1//2)/10,
                                  initstep = ODE.dtinit(F, y0, tspan, eps(T)^T(1//3)/10, eps(T)^T(1//2)/10, order=ODE.order(stepper))
                                  )
            solver = ODE.solve(equation,stepper,opts)
            solution = collect(ODE.dense(solver))
            nn = length(solution)
            t = Array(T,nn)
            y = Array(typeof(y0),nn)

            for (n,(tt,yy)) in enumerate(solution)
                t[n] = tt
                y[n] = isa(y0,Number) ? yy[1] : yy
            end
            y = hcat(y...).'
            @test maximum(abs(y-ana(t,y))) < tol
        end
    end
end

# Testing the lower-level iteration API

# Test sets with in-place F! [F!, y0, tspan, analytic]

function test_iterator_in_place()
    println("\nTesting iterators using in-place functions:")

    for rks in steppers
        println("Testing $rks")
        for ts in testsets_vector
            F!,y0,tspan,ana = ts
            T = eltype(tspan)
            stepper = rks{T}()
            # jac = ODE.forward_jacobian(F!,y0)
            # jac! = (t,y,J) -> copy!(J,jac(t,y))
            equation  = ODE.ExplicitODE(tspan[1],y0,F!)
            opts = ODE.Options{T}(
                                  tspan    = tspan,
                                  reltol   = eps(T)^T(1//3)/10,
                                  abstol   = eps(T)^T(1//2)/10,
                                  initstep = 0.001
                                  )
            solver = ODE.solve(equation,stepper,opts)
            solution = collect(ODE.dense(solver))
            nn = length(solution)
            t = Array(T,nn)
            y = Array(typeof(y0),nn)

            for (n,(tt,yy)) in enumerate(solution)
                t[n] = tt
                y[n] = isa(y0,Number) ? yy[1] : yy
            end
            y = hcat(y...).'
            @test maximum(abs(y-ana(t,y))) < tol
        end
    end
end


test_ode()
test_iterator_out_place()
test_iterator_in_place()
