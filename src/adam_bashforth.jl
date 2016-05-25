####################################
# Implicit Adam-Bashforth solvers
####################################
# (Hairer & Wanner 1996, Vol I, p.357-358


##Function: ode_imp_ab(F,y0, tspan, order(optional))
##order is set to 4 by default
##Adam Bashforth is a fixed step, fixed order, multistep Implicit method
function ode_imp_ab(F::Function,y0, tspan,order=4 ::Integer)
    if (0 <= order <= 3)
        coef = am_imp_coefficients3
        b = ms_coefficients4
    else
        #calculating higher order coefficients for implicit Adam Multon method
        coef = zeros(order+1, order+1)
        coef[1:4,1:4] = am_imp_coefficients3
        for s = 4:order
            for j = 0:s
                # Assign in correct order for multiplication below
                #  (a factor depending on j and s) .* (an integral of a polynomial with -(-1:s-1), except -(j-1), as roots)
                p_int = polyint(poly(diagm(-[-1:j - 2; j:s-1])))
                coef[s+1, s+1 - j] = ((-1)^j / factorial(j)
                               / factorial(s - j) * polyval(p_int, 1))
            end
        end
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

    h = diff(tspan)
    y = Array(typeof(y0), length(tspan))
    ydot = similar(y)
    y[1] = y0
    ydot[1] = F(tspan[1],y[1])

    ##PECP Method for Implicit Adam solver
    for i=1:length(tspan)-1
        steporder = min(i,order)

        ##(P)redict y[i+1] using explicit Adam Bashforth coefficients
        y[i+1]=y[i]
        for j=1:steporder
            y[i+1] += h[i]*b[steporder, j]*ydot[i-(steporder-1) + (j-1)]
        end

        ##(E)valuate function F at the approximate point tspan[i+1],y[i+1]
        ydot[i+1] = F(tspan[i+1],y[i+1])

        ##(C)orrect the formula using implicit Adam Multon coefficients
        y[i+1]=y[i]
        for j = 1 :steporder+1
            y[i+1] += h[i]*coef[steporder+1, j]*ydot[i - (steporder) + j]
        end

        ##(E)valulate the function anew at corrected approximatiion tspan[i+1],y[i+1]
        ydot[i+1] = F(tspan[i+1], y[i+1])
    end
    #println("Adam...Bash...nuff said")
    return tspan, y
end

##Adam Multon Coefficients for Implicit Method
const am_imp_coefficients3 =[1      0       0       0
                            1/2     1/2     0       0
                            -1/12    8/12    5/12   0
                            1/24    -5/24   19/24   9/24]
