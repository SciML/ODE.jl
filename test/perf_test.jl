# comparison of adaptive methods


# adaptive solvers
solvers = [
    ODE.ode23,
    ODE.ode23_bs,

    ODE.ode45_dp,
    ODE.ode45_fb,
    ODE.ode45_ck
    ]

println("\n== problem A3 in DETEST ==\n")
for solver in solvers
    println("testing method $(string(solver))")
    
    # problem A3 in DETEST
    t,y = solver((t,y)->y*cos(t), 1., [0.,20.]);
    @time begin
        for i=1:10
            t,y = solver((t,y)->y*cos(t), 1., [0.,20.]);
        end
    end

    println("*  number of steps: $(length(t))")
    println("*  minimal step: $(minimum(diff(t)))")
    println("*  maximal step: $(maximum(diff(t)))")
    
    println("*  maximal deviation from known solution: $(maximum(abs(y[:]-exp(sin(t)))))\n")
end

# problem B5 in DETEST
# motion of a rigid body without external forces
# (see also Example 1 from http://www.mathworks.se/help/matlab/ref/ode45.html)
function rigid(t,y)
	dy = copy(y)
	dy[1] = y[2] * y[3]
	dy[2] = -y[1] * y[3]
	dy[3] = -0.51 * y[1] * y[2]

	return dy
end

println("\n== problem B5 in DETEST ==\n")
for solver in solvers
    println("testing method $(string(solver))")
    
    # problem B5 in DETEST
    t,y = solver(rigid, [0., 1., 1.], [0.,12.]);
    @time begin
        for i=1:10
            t,y = solver(rigid, [0., 1., 1.], [0.,12.]);
        end
    end

    println("*  number of steps: $(length(t))")
    println("*  minimal step: $(minimum(diff(t)))")
    println("*  maximal step: $(maximum(diff(t)))\n")
    
    # println("*  maximal deviation from known solution: $(maximum(abs(y[:]-exp(sin(t)))))\n")
end
