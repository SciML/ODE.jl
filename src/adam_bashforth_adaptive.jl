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

    # initialization of step size
    tstart = tspan[1]
    tfinal = tspan[end]
    t[1] = tstart

    #step size array
    h = Float64[]
    h[1] = initstep
    if h[1] == 0.
        # initial guess at a step size
        h[1], tdir, F0 = hinit(F, y0, t, tfinal, 3, reltol, abstol)
    else
        tdir = sign(tfinal - t)
        F0 = F(t,y0)
    end
    #Note: h[n] = t[n+1]-t[n] => t[n+1] = t[n]+h[n]

    n = 1 #start at the second step (the first being y0)

    while abs(t - tfinal) > 0 && minstep < abs(h[n]) && successful_step = false
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
              c[j][q]=1/q
            elseif j==2
              c[j][q]=1/q/(q+1)
            else
              c[j][q] = c[j-1][q] - c[j-1][q+1]*h[n]/(t[n+1+1]-t[n-j+1+1])
            end
          end
        end

        #calculation of b_j(n), Φstar_j(n), and Φ_j(n) coefficients from recurence relations
        for j = 1:steporder
          if j== 1
            b[n][j] = 1
            Φstar[n][j] = dy[n]
            Φ[n][j] = dy[n]
            g[n][j] = c[j][1]
          else j>1
            #Note: since n>=j, when this runs, n>=2 (so n-1>=1) and
            b[n][j] = b[n][j-1]*(t[n+1+1]-t[n-j+1+1])/(t[n+1]-t[n-j+1])
            Φ[n][j+1] =  Φ[n][j]- Φstar[n-1][j]
            Φstar[n][j] = b[n][j]*Φ[n][j]
            g[n][j] = c[j][1]
          end
        end

        #Calculate next value: y[n+1].
        y[n+1] = y[n]
        for j= 1:steporder
          y[n+1] += h[n]*g[n][j]* Φstar[n][j]
        end

        ##(P)redict y[n+1]
        Φp[n+1][steporder+1] =
        g[steporder+1][]
        pred = y[n+1] + h[n]*g[n][steporder+1]*Φ[n+1][steporder+1]

        #error bounds dictating if attempted step was successful
        LE = h[n]*(g[n][steporder+1] - g[n][steporder])*Φp[n+1][steporder+1]

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
          h[1], tdir, F0 = hinit(F, y0, t, tfinal, 3, reltol, abstol)
      else
          tdir = sign(tfinal - t)
          F0 = F(t,y0)
      end
      #Note: h[n] = t[n+1]-t[n] => t[n+1] = t[n]+h[n]

      n = 1 #start at the second step (the first being y0)

      while abs(t - tfinal) > 0 && minstep < abs(h[n]) && successful_step = false
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
                c[j][q]=1/q
              elseif j==2
                c[j][q]=1/q/(q+1)
              else
                c[j][q] = c[j-1][q] - c[j-1][q+1]*h[n]/(t[n+1+1]-t[n-j+1+1])
              end
            end
          end

          #calculation of b_j(n), Φstar_j(n), and Φ_j(n) coefficients from recurence relations
          for j = 1:steporder
            if j== 1
              b[n][j] = 1
              Φstar[n][j] = dy[n]
              Φ[n][j] = dy[n]
              g[n][j] = c[j][1]
            else j>1
              #Note: since n>=j, when this runs, n>=2 (so n-1>=1) and
              b[n][j] = b[n][j-1]*(t[n+1+1]-t[n-j+1+1])/(t[n+1]-t[n-j+1])
              Φ[n][j+1] =  Φ[n][j]- Φstar[n-1][j]
              Φstar[n][j] = b[n][j]*Φ[n][j]
              g[n][j] = c[j][1]
            end
          end

          #Calculate next value: y[n+1] usin PECP Method
          ##(P)redict y[i+1] using explicit Adam Bashforth coefficients
          y[n+1]=y[n]
          for j=1:steporder
            #We use y[n+1] = y[n] + h[n] \sum_{j =0}^{k-1}(g[n][j]+ Φstar*[n][j])
            y[n+1] += h[n]*g[n][j]* Φstar[n][j]
          end
          pred = y[n+1] + h[n]*g[n][steporder+1]*Φ[n+1][steporder+1]

          ##(E)valuate function F at the approximate point t[i+1],y[i+1]
          dy[i+1] = F(t[i+1],pred)

          ##(C)orrect the formula using implicit Adam Multon coefficients
          #calculation of b_j(n), Φstar_j(n), and Φ_j(n) coefficients from recurence relations
          for j = 1:steporder
            if j== 1
              b[n][j] = 1
              Φstar[n][j] = dy[n]
              Φ[n][j] = dy[n]
              g[n][j] = c[j][1]
            else j>1
              #Note: since n>=j, when this runs, n>=2 (so n-1>=1) and
              b[n][j] = b[n][j-1]*(t[n+1+1]-t[n-j+1+1])/(t[n+1]-t[n-j+1])
              Φ[n][j+1] =  Φ[n][j]- Φstar[n-1][j]
              Φstar[n][j] = b[n][j]*Φ[n][j]
              g[n][j] = c[j][1]
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
          LE = h[n]*(g[n][steporder] - g[n][steporder-1])*Φp[n][steporder]

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
##Adam Multon Coefficients for Implicit Method
const am_imp_coefficients3 =[1      0       0       0
                            1/2     1/2     0       0
                            -1/12    8/12    5/12   0
                            1/24    -5/24   19/24   9/24]
