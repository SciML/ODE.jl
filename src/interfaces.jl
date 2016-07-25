"""

We assume that the initial data y0 is given at tspan[1], and that
tspan[end] is the last integration time.

"""

function ode{T,Y,S<:AbstractStepper}(F, y0::Y,
                                     tout::AbstractVector{T},
                                     stepper::Type{S};
                                     points = :all,
                                     kargs...)

    t0 = tout[1]

    # construct a solver
    equation  = explicit_ineff(t0,y0,F;kargs...)
    if points == :all
        solver = solve(equation, stepper;
                       tout = tout,
                       kargs...)
    elseif points == :specified
        solver = solve(equation, DenseStepper;
                       mehtod = stepper,
                       tout = tout,
                       kargs...)
    else
        error("Unsupported points value (should be :all or :specified)")
    end

    # determine if we have to unpack y
    extract = Y <: Number

    to = Array(T,0)
    yo = Array(Y,0)
    for (t,y) in solver
        push!(to,t)
        push!(yo, extract ? y[1] : copy(y))
    end

    return (to,yo)
end

"""
    ODE.odeXX(F,y0,t0;kargs...)

Solves an ODE `y'=F(t,y)` with initial conditions `y0` and `t0`.
"""

ode23s(F,y0,t0;kargs...)        = ode_conv(F,y0,t0,ModifiedRosenbrockStepper; kargs...)
ode1(F,y0,t0;kargs...)          = ode_conv(F,y0,t0,RKStepperFixed{:feuler}; kargs...)
ode2_midpoint(F,y0,t0;kargs...) = ode_conv(F,y0,t0,RKStepperFixed{:midpoint}; kargs...)
ode2_heun(F,y0,t0;kargs...)     = ode_conv(F,y0,t0,RKStepperFixed{:heun}; kargs...)
ode4(F,y0,t0;kargs...)          = ode_conv(F,y0,t0,RKStepperFixed{:rk4}; kargs...)
ode21(F,y0,t0;kargs...)         = ode_conv(F,y0,t0,RKStepperAdaptive{:rk21}; kargs...)
ode23(F,y0,t0;kargs...)         = ode_conv(F,y0,t0,RKStepperAdaptive{:rk23}; kargs...)
ode45_fe(F,y0,t0;kargs...)      = ode_conv(F,y0,t0,RKStepperAdaptive{:rk45}; kargs...)
ode45_dp(F,y0,t0;kargs...)      = ode_conv(F,y0,t0,RKStepperAdaptive{:dopri5}; kargs...)
const ode45 = ode45_dp
ode78(F,y0,t0;kargs...)         = ode_conv(F,y0,t0,RKStepperAdaptive{:feh78}; kargs...)


function ode_conv{Ty,T}(F,y0::Ty,t0::AbstractVector{T},stepper;kargs...)

    if !isleaftype(T)
        error("The output times have to be of a concrete type.")
    elseif !(T <:AbstractFloat)
        error("The time variable should be a floating point number.")
    end

    if !isleaftype(Ty) & !isleaftype(eltype(Ty))
        error("The initial data has to be of a concrete type (or an array)")
    end

    ode(F,y0,t0,stepper;kargs...)

end



"""

Convert a out-of-place explicitly defined ODE function to
ExplicitODE.  As the name suggests, the result is not going to be very
efficient.

"""
function explicit_ineff{T,Y}(t0::T, y0::AbstractVector{Y}, F::Function; kargs...)
    F!(t,y,dy) =copy!(dy,F(t,y))
    return ExplicitODE(t0,y0,F!; kargs...)
end

# A temporary solution for handling scalars, should be faster then the
# previous implementation.  Should be used only at the top level
# interface.  This function cheats by converting scalar functions F
# and jac to vector functions F! and jac!.  Still, solving this ODE
# will result in a vector of length one result, so additional external
# conversion is necessary.
function explicit_ineff{T,Y}(t0::T, y0::Y, F::Function; kargs...)
    F!(t,y,dy) =(dy[1]=F(t,y[1]))
    return ExplicitODE(t0,[y0],F!; kargs...)
end
