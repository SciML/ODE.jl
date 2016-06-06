####################################
# Explicit Adaptive Adam-Bashforth solvers
####################################
# (Hairer & Wanner 1996, Vol I, p.398-400, 421-423)
# ODE_MS Fixed-step, fixed-order multi-step numerical method
#   with Adams-Bashforth-Moulton coefficients
function ode_ms_adaptive(F, y0, tspan, order::Integer;reltol = 1.0e-5,
                                                      abstol = 1.0e-8,
                                                      minstep=abs(tspan[1] - tspan[1])/1e18,
                                                      maxstep=abs(tspan[end] - tspan[1])/2.5,
                                                      initstep=0.,
                                                      facmax = 1.5,
                                                      facmin = .1,
                                                      fac = .8)

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
    #test with fixed step size
    #@show h = diff(tspan)
    #Note: h[n] = t[n+1]-t[n] => t[n+1] = t[n]+h[n]

    #initialization of arrays used in calculation of further steps
    #all initialized to a 2x2 array, which will need to expand to nxn
    #over course of algorithm
    c = zeros(3,3)
    g = zeros(3,3)
    b = zeros(2,2)
    Φstar = Array{typeof(F(t[1], y[1])),2}(2,2)
    Φ = Array{typeof(F(t[1], y[1])),2}(2,2)

    LE = 1 #Error controller, defined here to be accessible throughout loop

    #start at n=1, the first step
    n = 1
    while abs(t[n] - tfinal) > 0 && sign(tfinal-t[n])==tdir #&&  minstep < abs(h[n])
        # Need to run the first several steps at reduced order
        if n == 1
            steporder = 1
        else
            steporder = min(n-1, order)
        end
        @show t[n]
        #Current derivative: set dy[n] = F(t[n], y[n])
        push!(dy, F(t[n], y[n]))

        #allocate space for t[n+1], y[n+1], dy[n+1], h[n+1]
        @show push!(t,t[1])
        push!(y,y[1])
        push!(dy,dy[1])
        push!(h,h[1])

        #resize coefficients arrays every power of two steporder
        if log2(n)%1 == 0 && n>1

          temp = zeros(2*n+4,2*n+4)
          #temp[1:n,1:n] = c[1:n,1:n]
          c = temp

          temp = zeros(2*n+4,2*n+4)
         #temp[1:n,1:n] = g[1:n,1:n]
          g = temp

          temp = zeros(2*n+4,2*n+4)
          #temp[1:n,1:n] = b[1:n,1:n]
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
          #@show n
          #@show trys = trys + 1
          #calculation (n+1)-th time position from n-th step size
          t[n+1]=t[n]+h[n]
          #####################################################################
          #Calculating explicit approximatiion to y[n+1]
          #####################################################################
          #calculation of c_{j,q} coefficients for j=0,...,k-1; q = 0,...,k-1-(j-1)
          for j = 1:steporder
            for q = 1:(steporder)-(j-1)
              if j ==1
                c[j,q]=1/q
              elseif j==2
                c[j,q]=1/q/(q+1)
              else
                c[j,q] = c[j-1,q] - c[j-1,q+1]*h[n]/(t[n+1]-t[n-(j-1)+1])
              end
            end
            g[n,j]=c[j,1]
          end
          #calculation of b_j(n), Φstar_j(n), and Φ_j(n) coefficients from recurence relations
          for j = 1:steporder
            if j== 1
              b[n,j] = 1
              Φ[n,j] = dy[n]
              Φstar[n,j] = dy[n]
            else j>1
              b[n,j] = b[n,j-1]*(t[n+1]-t[n-(j-1)+1])/(t[n]-t[n-(j-1)])
              Φ[n,j] =  Φ[n,j-1]- Φstar[n-1,j-1]
              Φstar[n,j] = b[n,j]*Φ[n,j]
            end
          end

          #Explicit approximation to y[n+1]
          ##Calculate next value: y[n+1]. Equation (5.5) in Hairer et al, Vol 1
          y[n+1] = y[n]
          for j= 1:steporder
            y[n+1] += h[n]*g[n,j]* Φstar[n,j]
          end
          #fixed step size code tested and gives correct result
          ######################################################################
          #Adaptive step selection
          ######################################################################
          #error bounds dictating if attempted step was successful
          #See Hairer et al equations (7.3 - 7.4)

          #Prediction of dy[n+1]
          dy[n+1] = F(t[n+1],y[n+1])

          #Don't use adaptive steps until n = order
          if n<order
            successful_step = true
            h[n+1]=tspan[n]
          else
          #Once we have enough initial points, we use adaptive step selection
          ##for this we need to calculate Φ[n+1,steporder+2]
            for j = steporder+1 #k
              if j== 1#0
                b[n,j] = 1
                Φ[n,j] = dy[n]
                Φstar[n,j] = dy[n]
              else j>1#0
                b[n,j] = b[n,j-1]*(t[n+1]-t[n-(j-1)+1])/(t[n]-t[n-(j-1)])
                Φ[n,j] =  Φ[n,j-1]- Φstar[n-1,j-1]
                Φstar[n,j] = b[n,j]*Φ[n,j]
              end
            end
            for j = 1:steporder
              if j== 1
                Φ[n+1,j] = dy[n+1]
              else j>1
                Φ[n+1,j] =  Φ[n+1,j-1]- Φstar[n,j-1]
              end
            end
            Φ[n+1,steporder+2] = Φ[n+1,steporder+1] -Φstar[n,steporder+1]

            ##Also we need to calculate g[n,steporder+1], g[n,steporder+2]
            for j = 1:steporder+2#k+1
              for q = 1:(steporder+2)-(j-1)#1:(k+1)-(j-1) #!!
                if j ==0
                  c[j,q]=1/q
                elseif j==1
                  c[j,q]=1/q/(q+1)
                else
                  c[j,q] = c[j-1,q] - c[j-1,q+1]*h[n]/(t[n+1]-t[n-(j-1)+1])
                end
              end
              g[n,j] = c[j,1]
            end

            #From this we can calculate the error
            @show LE = h[n]*(g[n,steporder+2] - g[n,steporder+1])*Φ[n+1,steporder+2]

            if norm(LE,Inf)>1  || maxstep < abs(h[n])
              successful_step = false
              h[n] = h[n]/2
            else
              successful_step = true
              @show h[n+1] = h[n]*(1/norm(LE,Inf))^(1/(steporder + 1)) #note: steporder = k + 1
            end # if statement for successful step
          end#while loop over successful step or not
        end

      #next step size selection
      #@show norm(LE,Inf)
     # @show h[n+1] = h[n]*(1/norm(LE,Inf))^(1/(steporder + 2)) #note: steporder = k + 1

      #determine optimal next order

      #move on to the next step
      n= n+1
    end #while loop over steps

    return vcat(t), y
end


####################################
# Implicit Adaptive Adam-Bashforth solvers
####################################
# (Hairer & Wanner 1996, Vol I, p.398-400, 421-423)

function ode_imp_adap_ab(F, y0, tspan, order::Integer;reltol = 1.0e-5,
                                                      abstol = 1.0e-8,
                                                      minstep=abs(tspan[end] - tspan[1])/1e18,
                                                      maxstep=abs(tspan[end] - tspan[1])/2.5,
                                                      initstep=0.,
                                                      facmax = 1.5,
                                                      facmin = .1,
                                                      fac = .8)

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
        h[1] = diff(tspan)[1]
    else
        tdir = sign(tfinal - tstart)
    end
    #test with fixed step size
    #@show h = diff(tspan)

    #Note: h[n] = t[n+1]-t[n] => t[n+1] = t[n]+h[n]

    #initialization of arrays used in calculation of further steps
    #all initialized to a 2x2 array, which will need to expand to nxn
    #over course of algorithm
    c = zeros(order+2,order+2)
    g = zeros(order+2,order+2)

    b = zeros(order+2,order+2)
    Φstar = Array{typeof(F(t[1], y[1])),2}(order+2,order+2)
    Φ = Array{typeof(F(t[1], y[1])),2}(order+2,order+2)

    LE = Array(Float64,1)#Error controller, defined here to be accessible throughout loop
    LE2 = Array(Float64,1)

    #start at n=2, the second step
    n = 1
    while abs(t[n] - tfinal) > 0 && sign(tfinal-t[n])==tdir #&&  minstep < abs(h[n])

        @show t[n]
        #Current derivative: set dy[n] = F(t[n], y[n])
        push!(dy, F(t[n], y[n]))

        #allocate space for t[n+1], y[n+1], dy[n+1], h[n+1]
        push!(t,t[1])
        push!(y,y[1])
        push!(dy,dy[1])
        push!(h,h[1])


        # Need to run the first several steps at reduced order
        if n == 1
          steporder = 1
        else
          steporder = min(n-1, order)
        end

        #resize coefficients arrays every power of two steporder
        if log2(n)%1 == 0 && n>1

          temp = zeros(2*n+1,2*steporder + 2)
          temp[1:n-1,1:steporder] = g[1:n-1,1:steporder]
          g = temp

          temp = zeros(2*steporder+2,2*steporder + 2)
          #temp[1:steporder+2,1:steporder+2] = c[1:steporder+1,1:steporder+1]
          c = temp

          temp = zeros(2*n+1,2*steporder + 2)
          #temp[1:n-1,1:steporder] = b[1:n-1,1:steporder]
          b = temp

          #Note, this is not necessary. Only need to store two most recent
          #Will fix this after solver is functional
          temp = [Φstar[1,1] for i=1:2*n+1, j=1:2*steporder+2]
          temp[1:n-1,1:steporder] = Φstar[1:n-1,1:steporder]
          Φstar = temp

          temp = [Φ[1,1] for i=1:2*n+1,j=1:2*steporder+2]
          temp[1:n-1,1:steporder] = Φ[1:n-1,1:steporder]
          Φ = temp
          println("Expanded!")
        end

        successful_step = false
        trys = 0
        while successful_step == false
          @show n
          @show trys = trys + 1
          #calculation (n+1)-th time position from n-th step size
          t[n+1]=t[n]+h[n]
          #####################################################################
          #Calculating explicit approximatiion to y[n+1]
          #####################################################################
          #calculation of c_{j,q} coefficients for j=0,...,k-1; q = 0,...,k-1-(j-1)
          for j = 1: steporder #j = 0:k-1
            for q = 1: steporder - (j-1)#1:k-(j-1)
              if j ==1#0
                c[j,q]=1/q
              elseif j==2#1
                c[j,q]=1/q/(q+1)
              else
                c[j,q] = c[j-1,q] - c[j-1,q+1]*h[n]/(t[n+1]-t[n-(j-1)+1])
              end
            end
            g[n,j] = c[j,1]
          end
          #calculation of b_j(n), Φstar_j(n), and Φ_j(n) coefficients from recurence relations
          for j = 1:steporder #0:k-1
            if j== 1#0
              b[n,j] = 1
              Φ[n,j] = dy[n]
              Φstar[n,j] = dy[n]
            else j>1#0
              #Note: since n>=j, when this runs, n>=2 (so n-1>=1) and
              b[n,j] = b[n,j-1]*(t[n+1]-t[n-(j-1)+1])/(t[n]-t[n-(j-1)])#!!! Need to fix this
              Φ[n,j] =  Φ[n,j-1]- Φstar[n-1,j-1]
              Φstar[n,j] = b[n,j]*Φ[n,j]
            end
          end
          ######################################################################
          #PECE Method
          ######################################################################
          #Explicit approximation to y[n+1]
          ##Calculate next value: y[n+1]. Equation (5.5) in Hairer et al, Vol 1
          #(P)rediction of y[n+1]
          y[n+1] = y[n]
          for j= 1:steporder#0:k-1
            y[n+1] += h[n]*g[n,j]* Φstar[n,j]
          end

          #(E)valation of dy[n+1] at t[n+1],y[n+1] <-prediction
          dy[n+1] = F(t[n+1],y[n+1])

          #add new coefficients  based on this prediction
          ##computations with recurrence relations
          for j = 1:steporder+1
            if j== 1#0
              Φ[n+1,j] = dy[n+1]
            else j>1#0
              #Note: since n>=j, when this runs, n>=2 (so n-1>=1) and
              Φ[n+1,j] =  Φ[n+1,j-1]- Φstar[n,j-1]
            end
          end
          ###calculation of c_{j,q} coefficients to calculate g_{k}(n)
          for j = 1:steporder+1#k
            for q = 1:(steporder+1)-(j-1)#1:(k+1)-(j-1) #!!
              if j ==0
                c[j,q]=1/q
              elseif j==1
                c[j,q]=1/q/(q+1)
              else
                c[j,q] = c[j-1,q] - c[j-1,q+1]*h[n]/(t[n+1]-t[n-(j-1)+1])
              end
            end
            g[n,j] = c[j,1]
          end

          #(C)orrection of y[n+1]
          y[n+1] = y[n+1] + h[n]*g[n,steporder+1]*Φ[n+1,steporder+1]

          #(E)valuate again
          dy[n+1] = F(t[n+1],y[n+1])

          #fixed step size code tested and gives correct result
          #... bugs adaptive step selection
          ######################################################################
          #Adaptive step selection
          ######################################################################
          #error bounds dictating if attempted step was successful
          #See Hairer et al equations (7.3 - 7.4)
          #Prediction of dy[n+1]
          dy[n+1] = F(t[n+1],y[n+1])

          #Don't use adaptive steps until n = order
          if n<order
            successful_step = true
            h[n+1]=tspan[n]
          else

            for j = steporder+1 #k
              if j== 1#0
                b[n,j] = 1
                Φ[n,j] = dy[n]
                Φstar[n,j] = dy[n]
              else j>1#0
                b[n,j] = b[n,j-1]*(t[n+1]-t[n-(j-1)+1])/(t[n]-t[n-(j-1)])
                Φ[n,j] =  Φ[n,j-1]- Φstar[n-1,j-1]
                Φstar[n,j] = b[n,j]*Φ[n,j]
              end
            end
            Φ[n+1,steporder+2] = Φ[n+1,steporder+1] -Φstar[n,steporder+1]

            for j = 1:steporder+2#k+1
              for q = 1:(steporder+2)-(j-1)#1:(k+1)-(j-1) #!!
                if j ==0
                  c[j,q]=1/q
                elseif j==1
                  c[j,q]=1/q/(q+1)
                else
                  c[j,q] = c[j-1,q] - c[j-1,q+1]*h[n]/(t[n+1]-t[n-(j-1)+1])
                end
              end
              g[n,j] = c[j,1]
            end
            #See Harier pg 422-423 for what is going on here
            ##Calculate truncation error LE in two different ways and update norm of error
            errors1=h[n]*(g[n,steporder+2] - g[n,steporder+1])*Φ[n+1,steporder+2]
            for error in errors1
                push!(LE,error)
            end
            ###(we test with two different norms)
            @show norm(LE)
            @show norm(LE,Inf)

            errors2 = h[n]* γ(steporder+2)*Φ[n+1,steporder+1]
            for error in errors2
                push!(LE2,error)
            end
            @show norm(LE2)
            @show norm(LE2,Inf)

            ##Continuing as Harier suggests,
            @show h[n]
            @show (h[n]*(1/norm(LE,Inf))^(1/(steporder+3))) #note: steporder = k + 1
            @show (h[n]*(1/norm(LE2,Inf))^(1/(steporder+3))) #note: steporder = k + 1

            @show hnewfac1= (1/norm(LE,Inf))^(1/(steporder)) #note: steporder = k + 1
            @show hnewfac2 = (1/norm(LE2,Inf))^(1/(steporder)) #note: steporder = k + 1

            ###See Harier pg 423, which heavily relies of Harier pg 168, for what follows
            if norm(LE,Inf)>1  || maxstep < abs(h[n])
              successful_step = false
              @show h[n] = max(h[n]*max(min(facmax,min(fac*hnewfac1, fac*hnewfac2)),facmin),minstep)
            else
              successful_step = true
              @show h[n+1] = max(h[n]*max(min(facmax,min(fac*hnewfac1, fac*hnewfac2)),facmin),minstep)
            end # if statement for successful step
          end#while loop over successful step or not
        end


      #determine optimal next order...will work on after step selection is working

      #move on to the next step
      n= n+1
    end #while loop over steps

    return vcat(t), y
end

###############################################################################
#Gamma function from Hairer et al pg 359
#-used in adaptive step selection (See Harier pg 423)
###############################################################################
function γ(j)
    if j == 0
        return 1
    elseif j > 0

    return polyval(((-1)^j / factorial(j))*polyint(poly(diagm(collect(1-(j-1):1)))),1)
    else
        error("order must be positive integer")
    end
end
