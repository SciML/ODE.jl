# symplectic integration methods

## velocity verlet with fixed step-size
## integrates the system
##     dx/dt = v(t)
##     dv/dt = a(t, x) = F(t, x)/m
## with initial conditions x_0 and v_0 using the velocity-verlet method.
##
## The function takes as a first argument a function, a(t,x),
## of time and position (the acceleration). The second argument
## specifies a list of times for which a solution is requested. The last
## two arguments are the initial conditions x_0 and v_0.
##
## The output, (t,x,v), consists of the solutions x(t), v(t) at times t=tspan.
##
function verlet_fixed(a, tspan, x_0, v_0)
	
	if length(x_0) != length(v_0)
		error("Intial data x_0 and v_0 must have equal length.")
	end
	
	x = Array(eltype(x_0), length(tspan), length(x_0))
	v = Array(eltype(v_0), length(tspan), length(v_0))
	
	x[1,:] = x_0'
	v[1,:] = v_0'
	
	atmp = a(tspan[1], x_0).'
	
    for i = 1:(length(tspan)-1)
		h = tspan[i+1] - tspan[i]
				
		v[i+1,:] = v[i,:] + atmp*h/2
		x[i+1,:] = x[i,:] + v[i+1,:]*h
		atmp = a(tspan[i]+h, x[i+1,:].').'
		v[i+1,:] = v[i+1,:] + atmp*h/2
	end
	
	return tspan, x, v
end 

## velocity verlet with adaptive step-size
## integrates the system
##     dx/dt = v(t)
##     dv/dt = a(t, x) = F(t, x)/m
## with initial conditions x_0 and v_0 using the velocity-verlet method.
## The step-size is adapted accoording to an error estimate based on
## comparing the results of two half steps and one full step.
## (for some background see: http://www.theorphys.science.ru.nl/people/fasolino/sub_java/pendula/computational-en.shtml)
##
## The function takes as a first argument a function, a(t,x),
## of time and position (the acceleration). The second argument
## specifies a list of times for which a solution is requested. The last
## two arguments are the initial conditions x_0 and v_0.
##
## The output, (t,x,v), consists of the solutions x(t), v(t) at the 
## intermediate integration times t including tspan[1] and tspan[end].
##
function verlet_hh2(a, tspan, x_0, v_0; atol = 1e-5, norm=Base.norm)
	if length(x_0) != length(v_0)
		error("Intial data x_0 and v_0 must have equal length.")
	end

	tc = tspan[1]
	tf = tspan[end]
	h = tf - tc
	v = v_0
	x = x_0
	
	tout = tc
	xout = x_0.'
	vout = v_0.'

	while tc < tf
		vo = v
		xo = x
		a_tmp = a(tc, xo)
		
		# try a full step
		v = vo + a_tmp*h/2
		x1 = xo + v*h
		v1 = v + a(tc+h, x)*h/2
		
		# try two half steps
		h = h/2
		
		v = vo + a_tmp*h/2
		x = xo + v*h
		a_tmp = a(tc+h, x)
		v =  v + a_tmp*h/2

		v =  v + a_tmp*h/2
		x =  x + v*h
		v =  v + a(tc+2*h, x)*h/2
		
		# error estimate
		err = (8/7)*norm(x1 - x)
		hmax = 2*h*abs(atol/err)^(1/4)
		
		if hmax < h
			# repeat step with smaller step size
			h = 0.9*2*hmax
			v = vo
			x = xo
		else
			# accept last step
			tc = tc + 2*h
			h = 2*h
			
			tout = [tout; tc]
			xout = [xout; x.']
			vout = [vout; v.']
		end		
	end
	return tout, xout, vout
end

# use adaptive method by default
verlet(a, tspan, x_0, v_0; atol = 1e-5, norm=Base.norm)= verlet_hh2(a, tspan, x_0, v_0; atol=atol, norm=norm)
