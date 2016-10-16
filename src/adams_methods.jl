using ProgressMeter
"""
    t,y = ode_ms(odefun, y0, tspan, order=4) with tspan = [t0:dt:tfinal]

Fixed-Equal Steps Adams-Bashforth multi-step solver for nonstiff problems

(Main ref: Hairer & Wanner 1996, Vol I, p.359-360, 372)

Note: If method is run at order k, the first step is taken at order 1,
and the order increases each step until it is order k on the kth step.
Thus results may be of lower accuracy than expected.
"""

function ode_ms(F, y0, tspan, order::Integer)
    dt = diff(tspan)

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

        # Check step size is equal or very near equal
        if !isapprox(dt[i],dt[1])
            error("Step sizes are not of constant length")
        end

        y[i+1] = y[i]
        for j = 1:steporder
            y[i+1] += dt[i]*b[steporder, j]*dy[i-(steporder-1) + (j-1)]
        end
    end
    return vcat(tspan), y
end

# Use order 4 by default
ode4ms(F, y0, tspan) = ODE.ode_ms(F, y0, tspan, 4)
ode5ms(F, y0, tspan) = ODE.ode_ms(F, y0, tspan, 5)


"""
    t,y = ode_am((odefun, y0, tspan, order=4) with tspan = [t0:dt:tfinal]

Fixed-Equal Step Adams-Moulton PECE solver for nonstiff problems

(Main ref: Hairer & Wanner 1996, Vol I, p.359-360, 372)

#Note: this function works for both scalar-like and array-like input

#Note: If method is run at order k, the first step is taken at order 1,
and the order increases each step until it is order k on the kth step.
Thus results may be of lower accuracy than expected.
"""

##Function: ode_am(F,y0, t, order(optional))
##order is set to 4 by default
function ode_am(F,y0, t,order::Integer = 4)
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

    ##Predict-Evalulate-Correct-Evaluate (PECE) Method for Implicit Adams solver
    for i=1:length(t)-1
        steporder = min(i,order)

        # Check step size is equal or very near equal
        if !isapprox(dt[i],dt[1])
            error("Step sizes are not of constant length")
        end

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

ode4am(F::Function,y0, t)= ode_am(F::Function,y0, t,4)
ode5am(F::Function,y0, t)= ode_am(F::Function,y0, t,5)

###############################################################################
# TABLEAUS FOR ADAMS FIXED STEP SOLVERS
###############################################################################

##Adams Bashforth Coefficients for Explicit Method
const ms_coefficients4 = [1        0        0         0
                          -1//2    3//2     0         0
                          5//12    -4//3    23//12    0
                          -9//24   37//24   -59//24   55//24]

##Adams Moulton Coefficients for Implicit Method
const am_imp_coefficients3 =   [1         0        0         0
                                1//2      1//2     0         0
                                -1//12    8//12    5//12     0
                                1//24     -5//24   19//24    9//24]


const γ_star_array=[1.0, -0.5,-0.08333333333333334,-0.041666666666666664,
                            -0.026388888888888885,-0.018750000000000003,
                            -0.014269179894179893,-0.01136739417989418,
                            -0.009356536596119928,-0.00789255401234568,
                            -0.006785849984634706,-0.005924056412337663,
                            -0.005236693257950285,-0.004677498407042264,
                            -0.004214952239005473,-0.003826899553211884,
                            -0.003497349845349918]
"""
#     t,y = ODE113(odefun, y0, tspan,order=4;kwords...)
#
#with tspan = [t0,tfinal] or[t0,t1,...,tfinal] for dense output
#
#Variable Step Variable Order (VSVO) Adams-Moulton Solver for nonstiff problems
#
#(Main ref: Hairer & Wanner 1996, Vol I, p.356-360, p.398-400, p.421-423)
"""
function ode113(F, y0, tspan; kwords...)
    # For y0 which don't support indexing.
    F_ = (t, y) -> [F(t, y[1])]
    t,y = ode113(F_, [y0], tspan; kwords...)
    return t, vcat_nosplat(y)
end

"a variable stepsize, variable order multi-step solver for nonstiff ODEs"
function ode113(F, y0::AbstractVector, tspan; order::Integer = 4,
                                        adaptive_order = true,
                                        max_order = 13,
                                        reltol = 1.0e-9,
                                        abstol = reltol,
                                        minstep = abs(tspan[1] - tspan[end])/1e18,
                                        maxstep = abs(tspan[end] - tspan[1])/2.5,
                                        initstep = 0.,
                                        fac = .8,
                                        norm = Base.norm,
                                        points = :all,
                                        progressmeter = false)

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
    # TODO: fix after https://github.com/JuliaLang/julia/issues/16667
    #y_n = deepcopy(y0)
    y_n = copy(y0)

    # initialzation of dy[n], which will have variable length as well
    dy_n = F(t[1],y0)
    Tdy = typeof(dy_n)

    # Declare work variables for t_np1, y_np1 and dy_np1  which are for
    # for t_{n+1}, y(t_{n+1}), and dy(t_{n+1})
    t_np1   = tspan[1]
    y_np1   = similar(y0, Eyf, dof)
    dy_np1  = similar(dy_n)

    # variables for dense output
    tout = [tstart]
    iter = 2
    yout = Array(Ty, 1)
    yout[1] = copy(y0)

    #initialization of step size dt, which will have variable length as well
    dt = initstep # current trial stepsize for a give step
    prev_dt = dt # the stepsize of the last successful step
    next_dt = 0.0 # stepsize to be used in next step

    if dt == 0.
        # initial guess at a step size
        # dt, tdir, dy0 = hinit(F, y0, tstart, tfinal, 3, reltol, abstol)
        dt, tdir, dy0 = hinit(F, y0, tstart, tfinal, 1, reltol, abstol)
        dt = dt
    else
        tdir = sign(tfinal - tstart)
    end

    n = 1     #`n` counts the time steps. We start at n=1, the first step
    num_stages = 0 # `num_stages`+1 is the order of the AB method at the current step
                   # Note: that due to Julia starting indexing at 1, whereas Hairer
                   # indices his variable for order, k, from 0, we have
                   # k= num_stages in what follows

    # initialization of arrays used in estimation y_np1 given dy[n],y[n], dt and
    # list order timestep are  initialized to a maxorder x maxorder dimension
    c,g,ϕstar_n, ϕ_n, ϕstar_nm1, b, ϕ_np1 = intialize_coef_arrays(max_order+3,Tdy,Et)

    if max_order <= 13
        γ_star = γ_star_array
    else
        γ_star= gamma_star_constr(max_order+3)
    end

    # Declaration of timeout variables
    timeout_VS_const = 0 # How many steps before a variable step size can be taken
    timeout_VS = timeout_VS_const # `timeout_VS` is a counter for this in the loop

    timeout_VO_const = 2 # How many steps before a variable step order can be taken
    timeout_VO = timeout_VO_const # `timeout_VO` is a counter for this in the loop

    consec_const_steps = 0 # counts how many consecutives steps have been taken
                           # with same step size

    # Declaration of  four local error (LE) or truncation errors necessary
    # for variable order calculations, which correspond to infinity norms of:
    #      LE(n+1)_{k-1}, LE(n+1)_{k}, E(n+1)_{k+1}, and LE(n+1)_{k+2}
    # from Hairer pgs 422-424
    err_km1, err_k, err_kp1, err_kp2 = 0.0,0.0,0.0,0.0

    #-##########################################################################
    # Variable Step Variable Order (VSVO) Step Integration Loop
    #-##########################################################################
    progressmeter
    if progressmeter
        prog = Progress(round(Int,1000*tdir*(tfinal-tstart)), "Solving ODE:")
    end
    while  tfinal*tdir > t[n]*tdir
        #--#####################################################################
        # Attempt a step with current trial stepsize `dt`
        #--#####################################################################
        successful_step = false
        while !successful_step
            #---################################################################
            # Make an approriate time step with current trial stepsize `dt`
            #---################################################################
            #Make the step
            t_np1 = t[n]+dt
            #---################################################################
            # Calculating explicit approximatiion to y_np1 with current step dt
            #---################################################################

            # Calculation of g_{j}(n)coefficients for j=0,...,k+1;
            g_coefs!(g,c,dt,t,t_np1,n,min(n,(num_stages+1)+1))
            #@show min(max(n-1,1),(num_stages)+1)
            # Calculation of  ϕstar_{j}(n), and ϕ_{j}(n) coef. for j=0,...,k;
            #ϕ_slots = min(max(n-1,1),(num_stages)+1)
            ϕ_and_ϕstar_coefs!(ϕstar_n, ϕ_n, ϕstar_nm1, b,dy_n,t,t_np1,n,1:(num_stages)+1)
            #Note that  1<= min(max(n-1,1),num_stages+1) <= n-1 and num_stages+1

            # Calculate predictor with explicit method of order `num_stages`. See
            # equation (5.5) in Hairer et al, Vol 1, where we note that here
            # y_np1 = y_{k}(t_{n+1})

            copy!(y_np1,y_n)
            for j= 0:num_stages-1
                y_np1 += dt*g[(j)+1]*ϕstar_n[(j)+1]
            end

            # Compute current estimate for derivative
            dy_np1=F(t_np1,y_np1)

            #---################################################################
            # Correction-Evaluaion (CE) Loop for implicit solver
            #     See Hairer p. 360 for a general description of this method
            #---################################################################

            # Set CE_loops to be how many correction-evaulation cycles
            # desired. As least one loop is necessary, for the sake of the
            # adaptive step selection
            CE_loops = 1
            for loop_counter=1:CE_loops
                # (C)orrection of y_np1
                # Add new coefficients  based on this prediction
                # computation of ϕ, ϕstar using recurrence relations
                ϕ_coefs!(ϕ_np1, ϕstar_n,dy_np1,((0)+1):((num_stages)+1))

                # Calculation of g_{num_stages + 1}(n) coefficients using
                # recurrence relations has already been done above, and does
                # not depend on dy_np1 or y_np1, so no update necessary

                # We correct the computation of y_np1 using y_np1
                # as part of our prediction Hairer equation 5.7
                y_np1 += dt*g[(num_stages)+1]*ϕ_np1[(num_stages)+1]

                # (E)valuate estimate for next derivative again
                dy_np1=F(t_np1,y_np1)
            end
            #---################################################################
            # Error Check for step size `dt` and Variable next step selection:
            #     Check error to determine if step size should be accepted:
            #     ->if yes, compute `next_dt` for next step
            #     ->if no, use this `next_dt`` for next attempt at this step
            #---################################################################

            #Don't allow increases of steps until n >= num_stages + 2
            if n==1
                successful_step = true
                next_dt = dt/4
            else
                # We calculate local truncation error, LE_{k+1}; the difference
                # between implicit results for order k and k+1.

                # We need ϕ_{k+1}(n+1) and g_{k+1}(n);the latter has already
                # been calculated above, the former still needs calculation
                #@show n, num_stages
                ϕ_coefs!(ϕ_np1, ϕstar_n,dy_np1,(num_stages+1)+1)

                # Calculate ||LE_{k+1}|| and next stepsize

                yerr_kp1= dt*(g[(num_stages+1)+1] - g[(num_stages)+1])*ϕ_np1[(num_stages+1)+1]
                err_kp1, next_dt, timeout_VS = stepsize_hw92!(dt,tdir,y_n,y_np1,
                                            yerr_kp1, num_stages, timeout_VS,dof,
                                            abstol, reltol, maxstep, norm)

                # Hairer equation 7.4: If ||LE_{k+1}|| =< 1
                # then step size is accepted and use next_dt for next step
                # Else, then step size is rejected, reattempt with fac*next_dt
                if err_kp1 <= 1
                    successful_step = true
                else
                    successful_step = false
                    dt = fac*next_dt
                end
            end # End of error check for current step size dt
        end # End of loop for finding successfull step size for nth step
            # repeats if success_step == false

        #Keep step size constant if `dt` is close to `prev_dt`
        if abs(next_dt/prev_dt)<2 && abs(next_dt/prev_dt)>=1
            consec_const_steps += 1
            next_dt = prev_dt
        else
            consec_const_steps = 0
        end

        #Don't step past `tfinal`
        if tdir*(t_np1 + next_dt) >= tdir*(tfinal)
            next_dt = tfinal-t_np1
        end
        push!(t,copy(t_np1))

        #--#####################################################################
        # Interpolation for dense output
        #--#####################################################################

        # interpolate specified points
        # TODO: less expensive interpolation. Core problem is that current
        # frame work does not allow for inexpensive interpolation for uneven
        # step sizes.

        while iter <= length(tspan) && tdir*tspan[iter]<tdir*t_np1
            dt_iter = tspan[iter]-t[n]

            # Calculation of g_{j}(n)coefficients for j=0,...,k+1;
            g_coefs!(g,c,dt,t,tspan[iter],n,min(n,num_stages+1))

            # Calculation of  ϕstar_{j}(n), and ϕ_{j}(n) coef. for j=0,...,k;
            ϕ_and_ϕstar_coefs!(ϕstar_n, ϕ_n, ϕstar_nm1, b,dy_n,t,tspan[iter],n,1:min(max(n-1,1),num_stages+1))

            yout_iter = copy(y_n)
            for j = 0:num_stages-1
                yout_iter += dt_iter*g[(j)+1]*ϕstar_n[(j)+1]
            end

            # CE iteration for more precise interpolation
            if n>=num_stages+1
                dy_iter=F(tspan[iter],yout_iter)
                ϕ_coefs!(ϕ_np1, ϕstar_n,dy_iter,1:(num_stages)+1)

                yout_iter += dt_iter*g[(num_stages)+1]*ϕ_np1[(num_stages)+1]
            end

            push!(tout,tspan[iter])
            push!(yout, copy(yout_iter))
            iter += 1
        end

        if iter <= length(tspan) && tspan[iter] == t_np1
            # handling when tspan and choosen time steps overlap
            push!(tout, copy(t_np1))
            push!(yout, copy(y_np1))
            iter += 1
        elseif points == :all
            push!(tout, copy(t_np1))
            push!(yout, copy(y_np1))
        end

        # restore values of coefficients for choosen step size
        g_coefs!(g,c,dt,t,t_np1,n,min(n,(num_stages)+1))
        ϕ_and_ϕstar_coefs!(ϕstar_n, ϕ_n, ϕstar_nm1, b, dy_n,t,t_np1,n,1:((num_stages)+1))
        ϕ_coefs!(ϕ_np1, ϕstar_n,dy_np1,1:((num_stages)+1))

        #--#####################################################################
        # Variable order selection for next step
        #   Not run until a success step was taken for n-th step
        #   See Hairer pg 424 for indepeth description
        #--#####################################################################

        current_stageorder = num_stages
        if n<=order-1
            # Need to run the first few steps at reduced order for it requires
            # the last 'order' number of points
            num_stages = min(n, max_order)
        elseif adaptive_order && num_stages>1
            # First, compute truncation errors for implicit orders: k-1,k
            # against current k+1

            yerr_km1 = dt*(g[(num_stages-1)+1] - g[(num_stages-2)+1])*ϕ_np1[(num_stages-1)+2]
            err_km1 = stepsize_hw92!(dt, tdir, y_n, y_np1, yerr_km1, num_stages,
                            timeout_VS,dof, abstol, reltol, maxstep, norm)[1]

            yerr_k = dt*(g[(num_stages)+1] - g[(num_stages-1)+1])*ϕ_np1[(num_stages)+1]
            err_k = stepsize_hw92!(dt, tdir, y_n, y_np1, yerr_k, num_stages,
                            timeout_VS, dof, abstol, reltol, maxstep, norm)[1]

            # Decrement num_stages if LE_{k} and LE_{k-1} < LE_{k+1}
            if max(err_km1,err_k) <= err_kp1
                num_stages = max(order,num_stages - 1)

            elseif consec_const_steps >= num_stages && timeout_VO==0 && n>(num_stages + 1)+1
                # Given enough consecutive steps taken at a constant step size
                # we compute truncation error for order implicit k+2
                ϕ_and_ϕstar_coefs!(ϕstar_n, ϕ_n, ϕstar_nm1, b, dy_n,t,t_np1,n,((num_stages+1)+1))
                ϕ_coefs!(ϕ_np1, ϕstar_n,dy_np1,(num_stages+2)+1)

                yerr_kp2 = dt*γ_star[(num_stages + 2)+1]*ϕ_np1[(num_stages+2)+1]
                err_kp2 = stepsize_hw92!(dt,tdir,y_n,y_np1,yerr_kp2,num_stages,
                            timeout_VS, dof, abstol, reltol, maxstep, norm)[1]

                # Increment num_stages if LE_{k+2} < LE_{k+1} and max(LE_{k},
                # LE_{k-1}) > LE_{k+1}
                if err_kp2 < err_kp1
                    num_stages = min(max_order,num_stages + 1)
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
        # Note: ϕ_n[:] and ϕ_star[2,:] corresponds to ϕ_n and ϕ*_n, and
        # ϕ_star[1,:] to ϕ*_{n-1}, etc

        ϕstar_nm1 = copy(ϕstar_n)
        #ϕstar_nm1[1:min(max(n-1,1),current_stageorder+1)] =
        #                        ϕstar_n[1:min(max(n-1,1),current_stageorder+1)]

        # Update stepsize variables
        prev_dt = dt
        dt = next_dt
        # push work variables
        y_n = y_np1
        dy_n = dy_np1

        # Increment our step
        n= n+1
        if progressmeter
          prog.counter = max(1,round(Int,1000.0*tdir*(t[n]-tstart)))
          ProgressMeter.updateProgress!(prog; showvalues = [(:steps,n),(:current_t,t[n]),(:tfinal,tfinal)])
        end
    end #while loop over steps

    return vcat(tout), yout
end

###############################################################################
#
# HELPER FUNCTIONS FOR ADAMS VSVO SOLVER
#
###############################################################################

###############################################################################
# Gamma star function from Hairer et al pg 359, table 1.2
# used in adaptive step selection (See Hairer pg 423)
###############################################################################
#precompute factorial_values
"compute γ*(j) for higher order versions of ode113"
function gamma_star_constr(index)
    γ_star = zeros(index)
    #γ_star[1] = 1
    for j = 1:index
        γ_star[j] = Polynomials.polyval((1 / factorial(j))*
            Polynomials.polyint(Polynomials.poly(diagm(collect(1-(j-1):1)))),1)
    end
    return γ_star
end
################################################################################
# Helper functions for calculating coefficients from recurence relationships
################################################################################
#Based on Hairer Lemma 5.1 on page 399
#Calculates c_j,q and g_n,j for n and j =1:num_stages
"compute c_{j,q} and g_{n,j} coefficients for ode113"
function g_coefs!(g,c,dt,t,t_np1,n,num_stages)
    for j = 1:num_stages
      for q = 1:((num_stages)-(j-1))
        if j ==1
          c[j,q]=1/q
        elseif j==2
          c[j,q]=1/q/(q+1)
        else
          c[j,q] = c[j-1,q] - c[j-1,q+1]*dt/(t_np1-t[n-(j-1)+1])
        end
      end
      g[j]=c[j,1]
    end
    return nothing
end


#Based on Hairer equation 5.9. Calculates b, ϕ, ϕstar for
#n and j  in range_of_indices
"compute ϕ_{n,j} and ϕ*_{n,j} coefficients for ode113"
function ϕ_and_ϕstar_coefs!(ϕstar_n, ϕ_n, ϕstar_nm1, b,dy_n,t,t_np1,n,range_of_indices)
    for j in range_of_indices #0:k-1
      if j== 1
        b[j] = 1
        ϕ_n[j] = copy(dy_n)
        ϕstar_n[j] = copy(dy_n)
      else
        b[j] = b[j-1]*(t_np1-t[n-(j-1)+1])/(t[n]-t[n-(j-1)])
        ϕ_n[j] =  ϕ_n[j-1]- ϕstar_nm1[j-1]
        ϕstar_n[j] = b[j]*ϕ_n[j]
      end
    end
    return nothing
end


#Based on Hairer equation 5.9. Just calculates the ϕ coefficients
"compute ϕ_{n+1,j} coefficients for ode113"
#ϕ_np1 as first
function ϕ_coefs!(ϕ_np1, ϕstar_n,dy_n,range_of_indices)
    for j  in range_of_indices
      if j== 1#0
          ϕ_np1[j] = copy(dy_n)
      else
        #Note: since n>=j, when this runs, n>=2 (so n-1>=1) and
        #@show ϕstar_n[j-1]
        ϕ_np1[j] = ϕ_np1[j-1] -  ϕstar_n[j-1]
      end
    end
    return nothing
end

################################################################################
# Helper functions initializing/resizing arrays
################################################################################
"used to initialize working arrays used in ode1113"
function intialize_coef_arrays(max_order,T,Et)
    c = zeros(Et,max_order,max_order)
    g = zeros(max_order)
    b = zeros(max_order)
    ϕstar_nm1 = Array{T}(1,max_order)
    ϕstar_n = Array{T}(1,max_order)
    ϕ_n = Array{T}(1,max_order)
    ϕ_np1 = Array{T}(1,max_order)
    return c,g,ϕstar_n, ϕ_n, ϕstar_nm1, b, ϕ_np1
end
