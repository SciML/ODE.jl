# comparison of adaptive methods

# problem B5 in DETEST
# motion of a rigid body without external forces
# (see also Example 1 from http://www.mathworks.se/help/matlab/ref/ode45.html)
function rigid(t,y)
	dy = similar(y)
	dy[1] = y[2] * y[3]
	dy[2] = -y[1] * y[3]
	dy[3] = -0.51 * y[1] * y[2]

	return dy
end

# selected problems from DETEST
problems = (["A3", true, (t,y)->y*cos(t), 1., 0., 20., (t)->exp(sin(t))],
            {"B5", false, rigid, [0., 1., 1.], 0., 12., nothing},)

# adaptive solvers
solvers = [
    ODE.ode23,
    ODE.ode23_bs,
    
    ODE.ode45_dp,
    ODE.ode45_fb,
    ODE.ode45_ck
    ]

for problem in problems

    (pname, hassol, F, y0, t0, tf, soly) = problem
    println("\n== problem $(pname) in DETEST ==\n")

    for solver in solvers
        println("testing method $(string(solver))")

        t,y = solver(F, y0, [t0, tf]);
        tau = @elapsed begin
            for i=1:10
                t,y = solver(F, y0, [t0, tf]);
            end
        end

        println("*  elapsed time: $(tau)")
        println("*  number of steps: $(length(t))")
        println("*  minimal step: $(minimum(diff(t)))")
        println("*  maximal step: $(maximum(diff(t)))")

        if hassol
            println("*  maximal deviation from known solution: $(maximum(abs(y[:]-soly(t))))\n")
        else
            println("\n")
        end
    end
end
