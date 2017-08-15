using Isaac
using Base.Test

println("Running tests for `Isaac.jl`")
@testset "Isaac" begin

include("testproblems.jl")
include("correctness.jl")
include("performance.jl")
include("test_dlanczos.jl")
include("test_saddle1.jl")
include("test_minim.jl")

end
