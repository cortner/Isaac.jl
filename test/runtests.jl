using CTKSolvers
using Base.Test

println("Running tests for `CTKSolvers.jl`")
@testset "StabilisedNewtonKrylov" begin

# include("testproblems.jl")
# include("correctness.jl")
# include("performance.jl")
# include("test_dlanczos.jl")
include("test_saddle1.jl")

end
