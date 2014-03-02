# comparison of adaptive methods

# adaptive solvers
solvers = [
    ODE.ode23,
    ODE.ode23_bs,

    ODE.ode45_dp,
    ODE.ode45_fb,
    ODE.ode45_ck
    ]

for solver in solvers
    println("testing method $(string(solver))")
    
    # problem A3 in DETEST
    t,y = solver((t,y)->y*cos(t), [0.,20.], [1.]);
    @time begin
        for i=1:10
            t,y = solver((t,y)->y*cos(t), [0.,20.], [1.]);
        end
    end

    println("*  number of steps: $(length(t))")
    println("*  minimal step: $(minimum(diff(t)))")
    println("*  maximal step: $(maximum(diff(t)))")
    
    println("*  maximal deviation from known solution: $(maximum(abs(y[:]-exp(sin(t)))))\n")
end