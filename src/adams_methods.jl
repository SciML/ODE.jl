################################################################################
# Fixed Step Adams-Bashforth multi-step solver for nonstiff problems
#
#     [T,X] = ode_ms(ODEFUN, Y0, TSPAN,ORDER=4) with TSPAN = [T0:H:TFINAL]
#
#     Note: If method is run at order k, the first step is taken at order 1,
#     and the order increases each step until it is order k on the kth step.
#     Thus results may be of lower accuracy than expected.
#     (Main ref: Hairer & Wanner 1996, Vol I, p.359-360, 372)
################################################################################

function ode_ms(F, y0, tspan, order::Integer)
    h = diff(tspan)
    y = Array(typeof(y0), length(tspan))
    y[1] = y0

    if 1 <= order <= 4
        b = ms_coefficients4
    else
        b = zeros(order, order)
        b[1:4, 1:4] = ms_coefficients4
        for s = 5:order
            for j = 0:(s - 1)
                # Assign in correct order for multiplication below
                #  (a factor depending on j and s) .* (an integral of a polynomial with -(0:s), except -j, as roots)
                p_int = polyint(poly(diagm(-[0:j - 1; j + 1:s - 1])))
                b[s, s - j] = ((-1)^j / factorial(j)
                               / factorial(s - 1 - j) * polyval(p_int, 1))
            end
        end
    end

    # TODO: use a better data structure here (should be an order-element circ buffer)
    dy = similar(y)
    for i = 1:length(tspan)-1
        # Need to run the first several steps at reduced order
        steporder = min(i, order)
        dy[i] = F(tspan[i], y[i])

        y[i+1] = y[i]
        for j = 1:steporder
            y[i+1] += h[i]*b[steporder, j]*dy[i-(steporder-1) + (j-1)]
        end
    end
    return vcat(tspan), y
end

# Use order 4 by default
ode4ms(F, y0, tspan) = ODE.ode_ms(F, y0, tspan, 4)
ode5ms(F, y0, tspan) = ODE.ode_ms(F, y0, tspan, 5)


################################################################################
#  Fixed Step Adams-Moulton PECE solver for nonstiff problems
#
#     [T,X] = ODE_AM(ODEFUN, Y0, TSPAN,ORDER=4) with TSPAN = [T0:H:TFINAL]
#
#     Note: If method is run at order k, the first step is taken at order 1,
#     and the order increases each step until it is order k on the kth step.
#     Thus results may be of lower accuracy than expected.
#     (Main ref: Hairer & Wanner 1996, Vol I, p.359-360, 372)
################################################################################


##Function: ode_am(F,y0, t, order(optional))
##order is set to 4 by default
function ode_am(F::Function,y0, t,order::Integer = 4)
    if (1 <= order <= 4)
        b_imp = am_imp_coefficients3
        b_exp = ms_coefficients4
    else
        #calculating higher order coefficients for implicit Adams Moulton method
        b_imp = zeros(order, order)
        b_imp[1:4,1:4] = am_imp_coefficients3
        k = order - 1 # For explicit method, order = k+1
        for s = 4:k
            for j = 0:s
                # Assign in correct order for multiplication below
                #  (a factor depending on j and s) .* (an integral of a polynomial with -(-1:s-1), except -(j-1), as roots)
                p_int = Polynomials.polyint(Polynomials.poly(diagm(-[-1:j - 2; j:s-1])))
                b_imp[s+1, s+1 - j] = ((-1)^j / factorial(j)
                               / factorial(s - j) * Polynomials.polyval(p_int, 1))
            end
        end
        b_exp = zeros(order, order)
        b_exp[1:4, 1:4] = ms_coefficients4
        k = order # For implicit method, order = k
        for s = 5:k
            for j = 0:(s - 1)
                # Assign in correct order for multiplication below
                #  (a factor depending on j and s) .* (an integral of a polynomial with -(0:s), except -j, as roots)
                p_int = Polynomials.polyint(Polynomials.poly(diagm(-[0:j - 1; j + 1:s - 1])))
                b_exp[s, s - j] = ((-1)^j / factorial(j)
                               / factorial(s - 1 - j) * Polynomials.polyval(p_int, 1))
            end
        end
    end

    dt = diff(t)
    y = Array(typeof(y0), length(t))
    dy = similar(y)
    y[1] = y0
    dy[1] = F(t[1],y[1])

    ##PECE Method for Implicit Adams solver
    for i=1:length(t)-1
        steporder = min(i,order)

        ##(P)redict y[i+1] using explicit Adams Bashforth coefficients
        y[i+1]=y[i]
        for j=1:steporder
            y[i+1] += dt[i]*b_exp[steporder, j]*dy[i-(steporder-1) + (j-1)]
        end

        ##(E)valuate function F at the approximate point t[i+1],y[i+1]
        dy[i+1] = F(t[i+1],y[i+1])

        ##(C)orrect the formula using implicit Adams Moulton coefficients
        y[i+1]=y[i]
        for j = 1 : steporder
            y[i+1] += dt[i]*b_imp[steporder, j]*dy[i - (steporder -1) + j]
        end

        ##(E)valulate the function anew at corrected approximatiion t[i+1],y[i+1]
        dy[i+1] = F(t[i+1], y[i+1])
    end
    return t, y
end

ode_am5(F::Function,y0, t)= ode_am(F::Function,y0, t,5)


################################################################################
# Variable Step Variable Order (VSVO) Adams-Moulton Solver for nonstiff problems
#
#     [T,X] = ODE113(ODEFUN, Y0, TSPAN,ORDER=4;kwords...) with TSPAN = [T0:TFINAL]
#     or[T0,T1,...,TFINAL] for dense output
#
#     (Main ref: Hairer & Wanner 1996, Vol I, p.356-360, p.398-400, p.421-423)
################################################################################
function ode113(F, y0, tspan, order::Integer=4; kwords...)
    # For y0 which don't support indexing.
    F_ = (t, y) -> [F(t, y[1])]
    t,y = ode113(F_, [y0], tspan, order; kwords...)
    return t, vcat_nosplat(y)
end

"a variable stepsize, variable order multi-step solver for nonstiff ODEs"
function ode113(F, y0::AbstractVector, tspan, order::Integer=4;adaptive_order = true,
                                                      max_order = 13,
                                                      reltol = 1.0e-9,
                                                      abstol = 1.0e-9,
                                                      minstep=abs(tspan[1] - tspan[end])/1e18,
                                                      maxstep=abs(tspan[end] - tspan[1])/2.5,
                                                      initstep=0.,
                                                      facmax = 1.5,
                                                      facmin = .1,
                                                      fac = .8,
                                                      norm=Base.norm,
                                                      points =:all)

    #-##########################################################################
    # Initializations
    #-##########################################################################
    Et, Eyf, Ty = make_consistent_types(F, y0, tspan)
    dof = length(y0)

    # initialization of time span, which will have variable length
    tspan = convert(Vector{Et}, tspan)
    tstart = tspan[1]
    tfinal = tspan[end]
    t = [tstart]

    # initialzation of y[n], which will have variable length as well
    y = Array(Ty, 1)
    allocate!(y, y0, dof)
    # TODO: fix after https://github.com/JuliaLang/julia/issues/16667
    # y[1] = deepcopy(y0)
    y[1] = copy(y0)

    # initialzation of dy[n], which will have variable length as well
    Tdy = typeof(F(t[1], y[1]))
    dy = Tdy[]
    push!(dy,F(t[1], y[1]))

    # preallocate work variables for t_np1, y_np1 and dy_np1  which are for
    # for t_{n+1}, y(t_{n+1}), and dy(t_{n+1})
    t_np1   = tspan[1]
    y_np1   = similar(y0, Eyf, dof)
    dy_np1  = similar(dy[1])

    # variables for dense output

    tout = [tstart]
    iter = 2
    yout = Array(Ty, 1)
    allocate!(yout, y0, dof)
    # TODO: see line 177
    yout[1] = copy(y0)

    #initialization of step size h, which will have variable length as well
    h = initstep # current trial stepsize for a give step
    prev_h = h # the stepsize of the last successful step
    next_h = 0.0 # stepsize to be used in next step

    if h == 0.
        # initial guess at a step size
        # h, tdir, dy0 = hinit(F, y0, tstart, tfinal, 3, reltol, abstol)
        h, tdir, dy0 = hinit(F, y0, tstart, tfinal, 1, reltol, abstol)

    else
        tdir = sign(tfinal - tstart)
    end

    n = 1     #`n` counts the time steps. We start at n=1, the first step
    steporder=1 # `steporder` is the order of the AB method at the current step
                # Note: that due to Julia starting indexing at 1, whereas Hairer
                # indices his variable for order, k, from 0, we have
                # k = steporder + 1 in what follows

    # initialization of arrays used in estimation y_np1 given dy[n],y[n], h and
    # list order timestep are  initialized to a maxorder x maxorder dimmension
    c,g,b,ϕstar, ϕ = intialize_coef_arrays(max_order+3,Tdy)

    timeout_VS_const = 0 # How many steps before a variable step size can be taken
    timeout_VS = timeout_VS_const # `timeout_VS` is a counter for this in the loop

    timeout_VO_const = 2 # How many steps before a variable step order can be taken
    timeout_VO = timeout_VO_const # `timeout_VO` is a counter for this in the loop

    consec_const_steps = 0 # counts how many consecutives steps have been taken
                           # with same step size

    # The four local error (LE) or truncation error necessary for variable
    # order calculations are listed below, which correspond to infinity norms of:
    #      LE(n+1)_{k-1}, LE(n+1)_{k}, E(n+1)_{k+1}, and LE(n+1)_{k+2}
    # from Harier pgs 422-424
    err_km1, err_k, err_kp1, err_kp2 = 0.0,0.0,0.0,0.0

    #-##########################################################################
    # Variable Step Variable Order (VSVO) Step Integration Loop
    #-##########################################################################
    while  tfinal*tdir > t[n]*tdir
        #--#####################################################################
        # Attempt a step with current trial stepsize `h`
        #--#####################################################################
        successful_step = false
        while !successful_step
            #---################################################################
            # Make an approriate time step with current trial stepsize `h`
            #---################################################################
            #Keep step size constant if `h` is close to `prev_h`
            if abs(h/prev_h)<2 && abs(h/prev_h)>=1
                consec_const_steps += 1
                h = prev_h
            else
                consec_const_steps = 0
            end

            #Don't step past `tfinal`
            if tdir*(t[n] + h) >= tdir*(tfinal)
                h = tfinal-t[n]
            end

            #Make the step
            t_np1 = t[n]+h
            #---################################################################
            # Calculating explicit approximatiion to y_np1 with current step h
            #---################################################################

            # Calculation of g_{j}(n)coefficients for j=0,...,k+1;
            g_coefs!(g,c,h,t,t_np1,n,min(n,steporder+2))

            # Calculation of  ϕstar_{j}(n), and ϕ_{j}(n) coef. for j=0,...,k;
            ϕ_and_ϕstar_coefs!(b, ϕ, ϕstar,dy[n],t,t_np1,n,1:min(max(n-1,1),steporder+1))
            #Note that  1<= min(max(n-1,1),steporder+1) <= n-1 and steporder+1

            # Calculate next value, y_np1, with order = k = steporder + 1. See
            # equation (5.5) in Hairer et al, Vol 1, where we note that here
            # y_np1 = y_{k}(t_{n+1})

            y_np1 = y[n]
            for j= 1:steporder
              y_np1 += h*g[j]* ϕstar[2,j]
            end
            # TODO: figure out why the following does not work
            # for j= 1:steporder
            #     for i=1:length(y_np1)
            #         y_np1[i] += h*g[j]*ϕstar[2,j][i]
            #     end
            # end

            # Compute current estimate for derivative
            dy_np1=F(t_np1,y_np1)
            #TODO: figure out why the following doesn't work
            #copy!(dy_np1,F(t_np1,y_np1))

            #---################################################################
            # Correction-Evaluaion (CE) Loop for implicit solver
            #     See Harier p. 360 for a general description of this method
            #---################################################################

            # Set CE_loops to be how many correction-evaulation cycles
            # desired. As least one loop is necessary, for the sake of the
            # adaptive step selection
            CE_loops = 1
            loop_counter = 0

            while loop_counter<CE_loops && n>=steporder+1
                loop_counter += 1
                # (C)orrection of y_np1
                # Add new coefficients  based on this prediction
                # computation of ϕ, ϕstar using recurrence relations
                ϕ_coefs!(ϕ, ϕstar,dy_np1,1:steporder+1)

                # Calculation of g_{steporder + 1}(n) coefficients using
                # recurrence relations has already been done above, and does
                # not depend on dy_np1 or y_np1, so no update necessary

                # We correct the computation of y_np1 using y_np1
                # as part of our prediction Harier equation 5.7

                y_np1 = y_np1 + h*g[steporder+1]*ϕ[3,steporder+1]
                # TODO: figure out why following does not work instead
                # for i=1:length(y_np1)
                #     y_np1[i] += h*g[steporder+1]*ϕstar[3,steporder+1][i]
                # end

                # (E)valuate estimate for next derivative again
                dy_np1=F(t_np1,y_np1)
            end
            #---################################################################
            # Error Check for step size `h` and Variable next step selection:
            #     Check error to determine if step size should be accepted:
            #     ->if yes, compute `next_h` for next step
            #     ->if no, use this `next_h`` for next attempt at this step
            #---################################################################

            #Don't use adaptive steps until n >= steporder + 2
            if n==1
                #TODO: adaptive step selection for n<steporder +2
                successful_step = true
                next_h = h/4
            elseif n<steporder+2
                yerr_k = h*(g[steporder+1] - g[steporder])*ϕ[3,steporder+1]
                err_k, next_h,timeout_VS = stepsize_hw92!(h, tdir, y[n], y_np1, yerr_k, steporder,
                                timeout_VS, dof, abstol, reltol, maxstep, norm)
                if err_k <= 1
                    successful_step = true
                else
                    successful_step = false
                    h = fac*next_h # make sure fac < 1
                end
            else
                # We calculate local truncation error, LE_{k+1}; the difference
                # between implicit results for order k and k+1.

                # We need ϕ_{k+1}(n+1) and g_{k+1}(n);the latter has already
                # been calculated above, the former still needs calculation
                ϕ_coefs!(ϕ, ϕstar,dy_np1,steporder+2)

                # Calculate ||LE_{k+1}|| and next stepsize
                yerr_kp1 = h*(g[steporder+2] - g[steporder+1])*ϕ[3,steporder+2]
                err_kp1, next_h, timeout_VS = stepsize_hw92!(h,tdir,y[n],y_np1,
                                            yerr_kp1, steporder, timeout_VS,dof,
                                            abstol, reltol, maxstep, norm)

                # Harier equation 7.4: If ||LE_{k+1}|| =< 1
                # then step size is accepted and use next_h for next step
                # Else, then step size is rejected, reattempt with fac*next_h
                if err_kp1 <= 1
                    successful_step = true
                else
                    successful_step = false
                    h = fac*next_h # make sure fac < 1
                end
            end # End of error check for current step size h
        end # End of loop for finding successfull step size for nth step
            # repeats if success_step == false



        #--#####################################################################
        # Interpolation for dense output
        #--#####################################################################
        # interpolate specified points
        while iter <= length(tspan) && tdir*tspan[iter]<tdir*t_np1
            h_iter = tspan[iter]-t[n]

            # Calculation of g_{j}(n)coefficients for j=0,...,k+1;
            g_coefs!(g,c,h,t,tspan[iter],n,min(n,steporder+1))

            # Calculation of  ϕstar_{j}(n), and ϕ_{j}(n) coef. for j=0,...,k;
            ϕ_and_ϕstar_coefs!(b, ϕ, ϕstar,dy[n],t,tspan[iter],n,1:min(max(n-1,1),steporder+1))

            yout_iter = y[n]
            for j= 1:steporder
                yout_iter += h_iter*g[j]*ϕstar[2,j]
            end

            #CE iteration for more precise interpolation
            if n>=steporder+1
                dy_iter=F(tspan[iter],yout_iter)
                ϕ_coefs!(ϕ, ϕstar,dy_iter,1:steporder+1)
                yout_iter = yout_iter + h_iter*g[steporder+1]*ϕ[3,steporder+1]
            end

            push!(tout,tspan[iter])
            push!(yout, yout_iter)
            iter += 1
        end

        if iter <= length(tspan) && tspan[iter] == t_np1
            # handling when tspan and choose time steps overlap
            push!(tout, t_np1)
            push!(yout, y_np1)
            iter += 1
        elseif points == :all
            push!(tout, t_np1)
            push!(yout, y_np1)
        end

        # restore values of coefficients for choosen step size
        g_coefs!(g,c,h,t,t_np1,n,min(n,steporder+1))
        ϕ_and_ϕstar_coefs!(b, ϕ, ϕstar,dy[n],t,t_np1,n,1:min(max(n-1,1),steporder+1))
        ϕ_coefs!(ϕ, ϕstar,dy_np1,1:steporder+1)

        # push work variables
        push!(y,y_np1)
        push!(dy,dy_np1)
        push!(t,t_np1)

        #--#####################################################################
        # Variable order selection for next step
        #   Note run until a success step was taken for n-th step
        #   See Harier pg 424 for indepeth description
        #--#####################################################################

        current_steporder = steporder
        if n<=order+1
            # Need to run the first few steps at reduced order for it requires
            # the last 'order' number of points
            steporder = n==1 ? 1 : n-1

        elseif adaptive_order && steporder>1
            # First, compute truncation errors for implicit orders: k-1,k
            # against current k+1
            yerr_km1 = h*(g[steporder] - g[steporder-1])*ϕ[3,steporder]
            err_km1 = stepsize_hw92!(h, tdir, y[n], y_np1, yerr_km1, steporder,
                            timeout_VS,dof, abstol, reltol, maxstep, norm)[1]

            yerr_k = h*(g[steporder+1] - g[steporder])*ϕ[3,steporder+1]
            err_k = stepsize_hw92!(h, tdir, y[n], y_np1, yerr_k, steporder,
                            timeout_VS, dof, abstol, reltol, maxstep, norm)[1]

            # Decrement steporder if LE_{k} and LE_{k-1} < LE_{k+1}
            if max(err_km1,err_k) <= err_kp1
                steporder = max(order,steporder - 1)

            elseif consec_const_steps >= steporder && timeout_VO==0 && n>steporder + 3
                # Given enough consecutive steps taken at a constant step size
                # we compute truncation error for order implicit k+2
                ϕ_and_ϕstar_coefs!(b, ϕ, ϕstar,dy[n],t,t_np1,n,steporder+2)
                ϕ_coefs!(ϕ, ϕstar,dy_np1,steporder+3)

                yerr_kp2 = h*γ_star(steporder + 3)*ϕ[3,steporder+3]
                err_kp2 = stepsize_hw92!(h,tdir,y[n],y_np1,yerr_kp2,steporder,
                            timeout_VS, dof, abstol, reltol, maxstep, norm)[1]

                # Increment steporder if LE_{k+2} < LE_{k+1} and max(LE_{k},
                # LE_{k-1}) > LE_{k+1}
                if err_kp2 < err_kp1
                    steporder = min(max_order,steporder + 1)
                    timeout_VO = timeout_VO_const # reset timout_VO
                end
            end
            # Keep timeout_VO 0, or decrement
            timeout_VO= timeout_VO==0? 0 : timeout_VO-1
        end # end of adaptive order selection

        #--#####################################################################
        # Prepare for next step
        #--#####################################################################

        # Shift data from this step, from the "current step" position to the
        # "previous step" position for ϕstar and ϕ
        # TODO: implement with DataStructures.CircularBuffer
        ϕstar[1,1:min(max(n-1,1),current_steporder+1)] =
                                ϕstar[2,1:min(max(n-1,1),current_steporder+1)]
        ϕ[1,1:min(max(n-1,1),current_steporder+1)] =
                                ϕ[2,1:min(max(n-1,1),current_steporder+1)]

        # Update stepsize pointers
        prev_h = h
        h = next_h
        # Increment our step
        n= n+1
    end #while loop over steps

    if points == :steps
        return vcat(t), y
    else
        return vcat(tout), yout
    end
end

###############################################################################
#
# HELPER FUNCTIONS and TABLEAUS FOR ADAMS SOLVERS
#
###############################################################################
##Adams Bashforth Coefficients for Explicit Method
const ms_coefficients4 = Float64[1        0        0         0
                                         -1//2    3//2     0         0
                                         5//12    -4//3    23//12    0
                                         -9//24   37//24   -59//24   55//24]

##Adams Moulton Coefficients for Implicit Method
const am_imp_coefficients3 = Float64[1         0        0         0
                                            1//2      1//2     0         0
                                            -1//12    8//12    5//12     0
                                            1//24     -5//24   19//24    9//24]



###############################################################################
# Gamma star function from Hairer et al pg 359, table 1.2
# used in adaptive step selection (See Harier pg 423)
###############################################################################
"compute γ*(j) for ode113"
function γ_star(j)
    if j == 0
        return 1
    elseif j > 0

    return Polynomials.polyval((1 / factorial(j))*
            Polynomials.polyint(Polynomials.poly(diagm(collect(1-(j-1):1)))),1)
    else
        error("order must be positive integer")
    end
end

################################################################################
# Helper functions for calculating coefficients from recurence relationships
################################################################################
#Based on Harier Lemma 5.1 on page 399
#Calculates c_j,q and g_n,j for n and j =1:steporder
"compute g_{n,j} coefficients for ode113"
function g_coefs!(g,c,h,t,t_np1,n,steporder)
    for j = 1:steporder
      for q = 1:(steporder)-(j-1)
        if j ==1
          c[j,q]=1/q
        elseif j==2
          c[j,q]=1/q/(q+1)
        else
          c[j,q] = c[j-1,q] - c[j-1,q+1]*h/(t_np1-t[n-(j-1)+1])
        end
      end
      g[j]=c[j,1]
    end
end

#Based on Harier equation 5.9. Calculates b, ϕ, ϕstar for
#n and j  in range_of_indices
"compute ϕ_{n,j} and ϕ*_{n,j} coefficients for ode113"
function ϕ_and_ϕstar_coefs!(b, ϕ, ϕstar,dy_n,t,t_np1,n,range_of_indices)
    for j in range_of_indices #0:k-1
      if j== 1#0
        b[j] = 1
        ϕ[2,j] = dy_n
        ϕstar[2,j] = dy_n
      else j>1#0
        #Note: since n>=j, when this runs, n>=2 (so n-1>=1) and
        b[j] = b[j-1]*(t_np1-t[n-(j-1)+1])/(t[n]-t[n-(j-1)])
        ϕ[2,j] =  ϕ[2,j-1]- ϕstar[1,j-1]
        ϕstar[2,j] = b[j]*ϕ[2,j]
      end
    end
end


#Based on Harier equation 5.9. Just calculates the ϕ coefficients
"compute ϕ_{n+1,j} coefficients for ode113"
function ϕ_coefs!(ϕ, ϕstar,dy_n,range_of_indices)
    for j  in range_of_indices
      if j== 1#0
        ϕ[3,j] = dy_n
      else j>1#0
        #Note: since n>=j, when this runs, n>=2 (so n-1>=1) and
        ϕ[3,j] =  ϕ[3,j-1]- ϕstar[2,j-1]
      end
    end
end

################################################################################
# Helper functions initializing/resizing arrays
################################################################################
"used to initialize working arrays used in ode1113"
function intialize_coef_arrays(dim,T)
    c = zeros(dim,dim)
    g = zeros(dim)
    b = zeros(dim)
    ϕstar = Array{T,2}(2,dim)
    ϕ = Array{T,2}(3,dim)
    return c,g,b,ϕstar, ϕ
end


function make_consistent_types(fn, y0, tspan)
    # There are a few types involved in a call to a ODE solver which
    # somehow need to be consistent:
    #
    # Et = eltype(tspan)
    # Ey = eltype(y0)
    # Ef = eltype(Tf)
    #
    # There are also the types of the containers, but they are not
    # needed as `similar` is used to make containers.
    # Tt = typeof(tspan)
    # Ty = typeof(y0)              # note, this can be a scalar
    # Tf = typeof(F(tspan(1),y0))  # note, this can be a scalar
    #
    # Returns
    # - Et: eltype of time, needs to be a real "continuous" type, at
    #       the moment a AbstractFloat
    # - Eyf: suitable eltype of y and f(t,y)
    #   --> both of these are set to typeof(y0[1]/(tspan[end]-tspan[1]))
    # - Ty: container type of y0

    # Needed interface:
    # On components: /, -
    # On container: eltype, promote_type
    # On time container: eltype

    Ty = typeof(y0)
    Eyf = typeof(y0[1]/(tspan[end]-tspan[1]))

    Et = eltype(tspan)
    @assert Et<:Real
    if !(Et<:AbstractFloat)
        Et = promote_type(Et, Float64)
    end

    # if all are Floats, make them the same
    if Et<:AbstractFloat &&  Eyf<:AbstractFloat
        Et = promote_type(Et, Eyf)
        Eyf = Et
    end

    !isleaftype(Et) && warn("The eltype(tspan) is not a concrete type!  Change type of tspan for better performance.")
    !isleaftype(Eyf) && warn("The eltype(y0/tspan[1]) is not a concrete type!  Change type of y0 and/or tspan for better performance.")

    return Et, Eyf, Ty
end
