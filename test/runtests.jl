using Isaac
using Base.Test

isCI = haskey(ENV, "CI")
notCI = !isCI

println("Running tests for `Isaac.jl`")
@testset "Isaac" begin

include("testproblems.jl")
include("correctness.jl")
if notCI
    include("performance.jl")
end 
include("test_dlanczos.jl")
include("test_saddle1.jl")
include("test_minim.jl")

end


# Two codes to look at for devising tests against NLopt
# https://github.com/JuliaDiffEq/DiffEqParamEstim.jl/blob/master/test/lorenz_true_test.jl
# https://github.com/JuliaDiffEq/DiffEqParamEstim.jl/blob/master/test/tests_on_odes/optim_test.jl#L16
