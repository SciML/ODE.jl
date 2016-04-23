# m3: This should be a generic function so methods can be added if
#     other types have something else than NaN for out-of-domain.
#     Although, better make it an option.
isoutofdomain = isnan

# m3: to be consistent: h->dt

#TODO make it a function on ExplicitODE and Options
function hinit{T}(F, y0, t0::T, reltol, abstol; tstop = T(Inf), order = 1)
    # Returns first step size
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
    y1 = similar(y0)
    for d = 1:length(y1)
        y1[d] = y0[d]+h0*f0[d]
    end
    f1 = F(t0 + h0, y1)
    # estimate second derivative
    d2 = norm(f1 - f0, Inf)/(tau*h0)
    if max(d1, d2) <= 10*eps(T)
        h1 = max(eps(T)^(1/3)/10, h0/10^3)
    else
        pow = -(2 + log10(max(d1, d2)))/(order+1)
        h1 = 10^pow
    end
    return min(100*h0, h1, abs(tstop-t0))
end
