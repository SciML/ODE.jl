Various basic Ordinary Differential Equation solvers implemented in Julia.

[![Build Status](https://travis-ci.org/JuliaLang/ODE.jl.png)](https://travis-ci.org/JuliaLang/ODE.jl) [![Coverage Status](https://img.shields.io/coveralls/JuliaLang/ODE.jl.svg)](https://coveralls.io/r/JuliaLang/ODE.jl)

Pull requests are always highly welcome to fix bugs, add solvers, or anything else!

# API discussions
There are currently discussions about how the Julian API for ODE solvers should look like, and the current documentation is more like a wishlist than a documentation. The API has changed considerably since the 0.1 release, so be carefull when you upgrade after the next version is released. 
# Current status of the project

The current release, v0.1, contains the basic functionality that was moved here when the package was originally moved here from Base. Although quite poorly tested, at least some of the functionality, especially the `ode45` solver, is quite reliable. However, that version is almost entirely undocumented, and will probably stay that way.

Since then, quite a lot has happened in the package, and the best way to use current ODE.jl is by checking out the latest master with `Pkg.checkout("ODE")`. By doing so, you get access to a new, better API -- but be careful; several breaking changes have been introduced since v0.1, and a few more might happen before v0.2. Therefore, the best way to learn the current API is to read the source. (The documentation in [http://github.com/JuliaLang/ODE.jl/master/blobs/doc/api.md](doc/api.md) is to be regarded as a wishlist, where some but not all of the features have been implemented as of yet).

Currently, `ODE` exports the following adaptive solvers:

* `ode23`: 2nd order adaptive solver with 3rd order error control, using the Bogacki–Shampine coefficients
* `ode45`: 4th order adaptive solver with 5th order error control, using the Dormand Prince coefficients. Fehlberg and Cash-Karp coefficients are also available.
* `ode78`: 7th order adaptive solver with 8th order error control, using the Fehlberg coefficients

all of which have the following basic API:

    tout, yout = odeXX(F, y0, tspan)
    
to solve the explicit ODE defined by dy/dt = F(t,y). A few other solvers are also exported, see the source code for details.

# Need something long-term reliable right now?

See [the Sundials.jl package](https://github.com/julialang/sundials.jl), which provides wrappers for the excellent Sundials ODE solver library.
