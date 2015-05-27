# hb1dae reformulates the hb1ode example as a differential-algebraic equation (DAE) problem.
# See Matlab help and http://radio.feld.cvut.cz/matlab/techdoc/math_anal/ch_8_od8.html#670396

using ODE, DASSL
M = diagm([1.,1,1])
M[3,3] = 0.

function hb1dae(t, y,ydot)
    res = zeros(length(y))
    hb1dae!(res, y,ydot)
    return res
end

function hb1dae!(res, y,ydot)
    res[1] = -0.04*y[1] + 1e4*y[2]*y[3] - ydot[1]
    res[2] = 0.04*y[1]  - 1e4*y[2]*y[3] - 3e7*y[2]^2 - ydot[2]
    res[3] = y[1] + y[2] + y[3] - 1
    return nothing
end
function Jhb1dae!(res, y, ydot, a)
    res[1,1] = -0.04 - a
    res[1,2] = 1e4*y[3]
    res[1,3] = 1e4*y[2]

    res[2,1] = 0.04
    res[2,2] = -1e4*y[3] - 3e7*2*y[2] -a;
    res[2,3] = -1e4*y[2]

    res[3,1] = 1.
    res[3,2] = 1.
    res[3,3] = 1.

    return nothing
end
function Jhb1dae(t, y, ydot, a)
    res = zeros(length(y),length(y))
    Jhb1dae!(res, y, ydot, a)
    return res
end


tspan = [0, 4e6]
y0 = [1., 0, 0]

t,ydassl = dasslSolve(hb1dae, y0, tspan, jacobian=Jhb1dae)
@time t,ydassl = dasslSolve(hb1dae, y0, tspan, jacobian=Jhb1dae)
y = hcat(ydassl...);

## ROSW
t,yrosw = ode_rosw(hb1dae!, Jhb1dae!, y0, tspan)
@time t,yrosw = ode_rosw(hb1dae!, Jhb1dae!, y0, tspan)
yr = hcat(yrosw...);

println("Relative difference between DASSL vs ROSW:")
println(abs(ydassl[end]-yrosw[end])./abs(ydassl[end]))


# using Winston
# plot(t, y[1,:], xlog=true)
# oplot(t, y[2,:]*1e4, xlog=true)
# oplot(t, y[3,:], xlog=true)
# oplot(t, y[1,:], xlog=true)
