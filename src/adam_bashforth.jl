####################################
# Implicit Adam-Bashforth solvers
####################################
# (Hairer & Wanner 1996, Vol I, p.357-358


##Function: ode_imp_ab(F,y0, t, order(optional))
##order is set to 4 by default
##Adam Bashforth is a fixed step, fixed order, multistep Implicit method
function ode_imp_ab(F::Function,y0, t,order::Integer = 4)
    if (1 <= order <= 4)
        b_imp = am_imp_coefficients3
        b_exp = ms_coefficients4
    else
        #calculating higher order coefficients for implicit Adam Multon method
        b_imp = zeros(order, order)
        b_imp[1:4,1:4] = am_imp_coefficients3
        k = order - 1 # For explicit method, order = k+1
        for s = 4:k
            for j = 0:s
                # Assign in correct order for multiplication below
                #  (a factor depending on j and s) .* (an integral of a polynomial with -(-1:s-1), except -(j-1), as roots)
                p_int = polyint(poly(diagm(-[-1:j - 2; j:s-1])))
                b_imp[s+1, s+1 - j] = ((-1)^j / factorial(j)
                               / factorial(s - j) * polyval(p_int, 1))
            end
        end
        b_exp = zeros(order, order)
        b_exp[1:4, 1:4] = ms_coefficients4
        k = order # For implicit method, order = k
        for s = 5:k
            for j = 0:(s - 1)
                # Assign in correct order for multiplication below
                #  (a factor depending on j and s) .* (an integral of a polynomial with -(0:s), except -j, as roots)
                p_int = polyint(poly(diagm(-[0:j - 1; j + 1:s - 1])))
                b_exp[s, s - j] = ((-1)^j / factorial(j)
                               / factorial(s - 1 - j) * polyval(p_int, 1))
            end
        end
    end

    dt = diff(t)
    y = Array(typeof(y0), length(t))
    dy = similar(y)
    y[1] = y0
    dy[1] = F(t[1],y[1])

    ##PECP Method for Implicit Adam solver
    for i=1:length(t)-1
        steporder = min(i,order)

        ##(P)redict y[i+1] using explicit Adam Bashforth coefficients
        y[i+1]=y[i]
        for j=1:steporder
            y[i+1] += dt[i]*b_exp[steporder, j]*dy[i-(steporder-1) + (j-1)]
        end

        ##(E)valuate function F at the approximate point t[i+1],y[i+1]
        dy[i+1] = F(t[i+1],y[i+1])

        ##(C)orrect the formula using implicit Adam Multon coefficients
        y[i+1]=y[i]
        for j = 1 : steporder
            y[i+1] += dt[i]*b_imp[steporder, j]*dy[i - (steporder -1) + j]
        end

        ##(E)valulate the function anew at corrected approximatiion t[i+1],y[i+1]
        dy[i+1] = F(t[i+1], y[i+1])
    end
    return t, y
end

##Adam Multon Coefficients for Implicit Method
const am_imp_coefficients3 =[1      0       0       0
                            1/2     1/2     0       0
                            -1/12    8/12    5/12   0
                            1/24    -5/24   19/24   9/24]
