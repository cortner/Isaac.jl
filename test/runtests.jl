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
