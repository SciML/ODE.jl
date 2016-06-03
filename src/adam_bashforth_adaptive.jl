####################################
# Explicit Adaptive Adam-Bashforth solvers
####################################
# ODE_MS Fixed-step, fixed-order multi-step numerical method
#   with Adams-Bashforth-Moulton coefficients
function ode_ms_adaptive(F, y0, tspan, order::Integer;reltol = 1.0e-5,
                                                      abstol = 1.0e-8,
                                                      minstep=abs(tspan[end] - tspan[1])/1e18,
                                                      maxstep=abs(tspan[end] - tspan[1])/2.5,
                                                      initstep=0.)

    # initialization of time span, which will have variable length
    tstart = tspan[1]
    tfinal = tspan[end]
    t = typeof(tspan[1])[]
    push!(t,tstart)

    # initialzation of y[n] and dy[n]
    y = typeof(y0)[]
    push!(y,y0)
    dy = typeof(F(t[1], y[1]))[]

    #initialization of step size array, while will have variable length as well
    h = collect(initstep)
    if h[1] == 0.
        # initial guess at a step size
        h[1], tdir, dy0 = hinit(F, y0, tstart, tfinal, 3, reltol, abstol)
    else
        tdir = sign(tfinal - tstart)
    end
    #Note: h[n] = t[n+1]-t[n] => t[n+1] = t[n]+h[n]

    #initialization of arrays used in calculation of further steps
    #all initialized to a 2x2 array, which will need to expand to nxn
    #over course of algorithm
    c = zeros(2,2)
    g = zeros(2,2)
    b = zeros(2,2)
    Φstar = Array{typeof(F(t[1], y[1])),2}(2,2)
    Φ = Array{typeof(F(t[1], y[1])),2}(2,2)

    LE = 1 #Error controller, defined here to be accessible throughout loop

    #start at n=1, the first step
    n = 1
    while abs(t[n] - tfinal) > 0 && sign(tfinal-t[n])==tdir &&  minstep < abs(h[n])
        # Need to run the first several steps at reduced order
        steporder = min(n, order)
        @show t[n]
        #Current derivative: set dy[n] = F(t[n], y[n])
        push!(dy, F(t[n], y[n]))

        #allocate space for t[n+1], y[n+1], dy[n+1], h[n+1]
        push!(t,t[1])
        push!(y,y[1])
        push!(dy,dy[1])
        push!(h,h[1])

        #resize coefficients arrays every power of two steporder
        if log2(n)%1 == 0 && n>1

          temp = zeros(2*n+4,2*n+4)
          temp[1:n,1:n] = c[1:n,1:n]
          c = temp

          temp = zeros(2*n+4,2*n+4)
          temp[1:n,1:n] = g[1:n,1:n]
          g = temp

          temp = zeros(2*n+4,2*n+4)
          temp[1:n,1:n] = b[1:n,1:n]
          b = temp

          temp = [Φstar[1,1] for i=1:2*n+4, j=1:2*n+4 ]
          temp[1:n-1,1:n-1] = Φstar[1:n-1,1:n-1]
          Φstar = temp

          temp = [Φ[1,1] for i=1:2*n+4, j=1:2*n+4 ]
          temp[1:n-1,1:n-1] = Φ[1:n-1,1:n-1]
          Φ = temp
          println("Expanded!")
        end

        successful_step = false
        trys = 0
        while successful_step == false
          trys = trys + 1
          #calculation (n+1)-th time position from n-th step size
          t[n+1]=t[n]+h[n]

          #calculation of c_{j,q} coefficients
          for j = 1:steporder
            for q = 1:(steporder+1)+1-j
              if j ==1
                c[j,q]=1/q
              elseif j==2
                c[j,q]=1/q/(q+1)
              else
                c[j,q] = c[j-1,q] - c[j-1,q+1]*h[n]/(t[n+1]-t[n-j+1])
              end
            end
          end

          #calculation of b_j(n), Φstar_j(n), and Φ_j(n) coefficients from recurence relations
          for j = 1:steporder
            if j== 1
              b[n,j] = 1
              Φ[n,j] = dy[n]
              Φstar[n,j] = dy[n]

              g[n,j] = c[j,1]
            else j>1
              #Note: since n>=j, when this runs, n>=2 (so n-1>=1) and
              b[n,j] = b[n,j-1]*(t[n+1]-t[n-(j-1)+1])/(t[n]-t[n-(j-1)])#!!! Need to fix this
              #maybe instead  b[n,j] = b[n,j-1]*(t[n+1]-t[)n-j+1])/(t[n]-t[n-(j-1)])#!!! with j = 1: steporder + 1

              Φ[n,j] =  Φ[n,j-1]- Φstar[n-1,j-1]
              Φstar[n,j] = b[n,j]*Φ[n,j]

              g[n,j] = c[j,1]
            end
          end

          #Calculate next value: y[n+1].
          y[n+1] = y[n]
          for j= 1:steporder
            y[n+1] += h[n]*g[n,j]* Φstar[n,j]
          end

          ##(P)redict y[n+1]
          ###Prediction of y'
          dy[n+1] = F(t[n+1],y[n+1])

          ###add new coefficients  based on this prediction

          ####calculation of c_{j,q} coefficients
          for j = steporder:steporder+1
            for q = 1:(steporder+1)-j
              if j ==1
                c[j,q]=1/q
              elseif j==2
                c[j,q]=1/q/(q+1)
              else
                c[j,q] = c[j-1,q] - c[j-1,q+1]*h[n]/(t[n+1]-t[n-j+1])
              end
            end
          end
          g[n,steporder + 1] = c[steporder+1,1]


          ####computations with recurrence relations
          for j = 1:steporder+1
            if j== 1
              Φ[n+1,j] = dy[n+1]
            else j>1
              #Note: since n>=j, when this runs, n>=2 (so n-1>=1) and
              Φ[n+1,j] =  Φ[n+1,j-1]- Φstar[n,j-1]
            end
          end

          pred = y[n+1] + h[n]*g[n,steporder+1]*Φ[n+1,steporder+1]

          #error bounds dictating if attempted step was successful
          @show h[n]
          @show g[n,steporder+1]
          @show g[n,steporder]
          @show Φ[n+1,steporder+1]
          LE = h[n]*(g[n,steporder+1] - g[n,steporder])*Φ[n+1,steporder+1]

          if norm(LE,Inf)>1  || maxstep < abs(h[n])
            successful_step = false
            h[n] = h[n]/2
          else
            successful_step = true
          end # if statement for successful step
        end #while loop over successful step or not

      #next step size selection
      @show norm(LE,Inf)
      h[n+1] = h[n]*(1/norm(LE,Inf))^(1/steporder + 1) #note: k +1 = steporder

      #determine optimal next order

      #move on to the next step
      n= n+1
    end #while loop over steps

    return vcat(t), y
end

#=
####################################
# Implicit Adaptive Adam-Bashforth solvers
####################################
# (Hairer & Wanner 1996, Vol I, p.359-360, 372)


##Function: ode_imp_ab(F,y0, t, order(optional))
##order is set to 4 by default
##Adam Bashforth is a fixed step, fixed order, multistep Implicit method
function ode_imp_adap_ab(F, y0, tspan, order::Integer = 4;reltol = 1.0e-5,
                                                        abstol = 1.0e-8,
                                                        minstep=abs(tspan[end] - tspan[1])/1e18,
                                                        maxstep=abs(tspan[end] - tspan[1])/2.5,
                                                        initstep=0.)

      # initialization of step size
      tstart = tspan[1]
      tfinal = tspan[end]
      t[1] = tstart

      #step size array
      h = Float64[]
      h[1] = initstep
      if h[1] == 0.
          # initial guess at a step size
          h[1], tdir, dy[1] = hinit(F, y0, tstart, tfinal, 3, reltol, abstol)
      else
          tdir = sign(tfinal - t)
          dy[1] = F(t,y0)
      end
      #Note: h[n] = t[n+1]-t[n] => t[n+1] = t[n]+h[n]

      n = 1 #start at the second step (the first being y0)

      while abs(t - tfinal) > 0 && minstep < abs(h[n]) && successful_step == false
          # Need to run the first several steps at reduced order
          steporder = min(n, order)

          #Current derivative
          dy[n] = F(t[n], y[n])

          #calculation n-th step size
          #h[n]
          t[n+1]=t[n]+h[n]

          #calculation of c_{j,q} coefficients
          for j = 1:steporder
            for q = 1:steporder-j
              if j ==1
                c[j,q]=1/q
              elseif j==2
                c[j,q]=1/q/(q+1)
              else
                c[j,q] = c[j-1,q] - c[j-1,q+1]*h[n]/(t[n+1+1]-t[n-j+1+1])
              end
            end
          end

          #calculation of b_j(n), Φstar_j(n), and Φ_j(n) coefficients from recurence relations
          for j = 1:steporder
            if j== 1
              b[n,j] = 1
              Φstar[n,j] = dy[n]
              Φ[n,j] = dy[n]
              g[n,j] = c[j,1]
            else j>1
              #Note: since n>=j, when this runs, n>=2 (so n-1>=1) and
              b[n,j] = b[n,j-1]*(t[n+1+1]-t[n-j+1+1])/(t[n+1]-t[n-j+1])
              Φ[n,j+1] =  Φ[n,j]- Φstar[n-1,j]
              Φstar[n,j] = b[n,j]*Φ[n,j]
              g[n,j] = c[j,1]
            end
          end

          #Calculate next value: y[n+1] usin PECP Method
          ##(P)redict y[i+1] using explicit Adam Bashforth coefficients
          y[n+1]=y[n]
          for j=1:steporder
            #We use y[n+1] = y[n] + h[n] \sum_{j =0}^{k-1}(g[n,j]+ Φstar*[n,j])
            y[n+1] += h[n]*g[n,j]* Φstar[n,j]
          end
          pred = y[n+1] + h[n]*g[n,steporder+1]*Φ[n+1,steporder+1]

          ##(E)valuate function F at the approximate point t[i+1],y[i+1]
          dy[i+1] = F(t[i+1],pred)

          ##(C)orrect the formula using implicit Adam Multon coefficients
          #calculation of b_j(n), Φstar_j(n), and Φ_j(n) coefficients from recurence relations
          for j = 1:steporder
            if j== 1
              b[n,j] = 1
              Φstar[n,j] = dy[n]
              Φ[n,j] = dy[n]
              g[n,j] = c[j,1]
            else j>1
              #Note: since n>=j, when this runs, n>=2 (so n-1>=1) and
              b[n,j] = b[n,j-1]*(t[n+1+1]-t[n-j+1+1])/(t[n+1]-t[n-j+1])
              Φ[n,j+1] =  Φ[n,j]- Φstar[n-1,j]
              Φstar[n,j] = b[n,j]*Φ[n,j]
              g[n,j] = c[j,1]
            end
          end

          y[i+1]=y[i]
          for j = 1 : steporder
              y[i+1] += dt[i]*b_imp[steporder, j]*dy[i - (steporder -1) + j]
          end

          ##(E)valulate the function anew at corrected approximatiion t[i+1],y[i+1]
          dy[i+1] = F(t[i+1], y[i+1])
        end

          #error bounds dictating if attempted step was successful
          LE = h[n]*(g[n,steporder] - g[n,steporder-1])*Φp[n,steporder]

          if norm(LE,Inf)>1
            successful_step = false
            h[n] = h[n]/2
          else
            successful_step = true

            #next step size selection
            h[n+1] = h[n]*(1/norm(LE,Inf))^(1/steporder + 1) #note: k +1 = steporder

            #determine optimal next order

            #move on to the next step
            n= n+1

          end # if statement for successful step
      end #while loop over steps

      return vcat(t), y
end
=#
