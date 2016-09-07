using ODE
using ODE.ODETests
using Base.Test

@testset "ODE tests" begin
    include("iterators.jl")
    include("top-interface.jl")
end
# TODO: do we still need this?
# include("interface-tests.jl")
