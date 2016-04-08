isoutofdomain = isnan

function hinit{T}(F, y0, t0::T, reltol, abstol; tstop = Inf, order = 1)
    # Returns first step size
    tdir = sign(tstop-t0)
    tau = max(reltol*norm(y0, Inf), abstol)
    d0 = norm(y0, Inf)/tau
    f0 = F(t0, y0)
    d1 = norm(f0, Inf)/tau
    if min(d0,d1) < eps(T)^(1/3)
        h0 = eps(T)^(1/3)/10
    else
        h0 = (d0/d1)/100
    end
    # perform Euler step
    y1 = y0 + tdir*h0*f0
    f1 = F(t0 + tdir*h0, y1)
    # estimate second derivative
    d2 = norm(f1 - f0, Inf)/(tau*h0)
    if max(d1, d2) <= 10*eps(T)
        h1 = max(eps(T)^(1/3)/10, h0/10^3)
    else
        pow = -(2 + log10(max(d1, d2)))/(order+1)
        h1 = 10^pow
    end
    return min(100*h0, h1, tdir*abs(tstop-t0))
end
