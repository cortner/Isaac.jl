using Isaac
using Test
using LinearAlgebra, SparseArrays, Random

print_tf(::Test.Pass) = printstyled("+", bold=true, color=:green)
print_tf(::Test.Fail) = printstyled("-", bold=true, color=:red)
print_tf(::Tuple{Test.Error,Bool}) = printstyled("x", bold=true, color=:magenta)


Random.seed!(1)

isCI = haskey(ENV, "CI")
notCI = !isCI

println("Running tests for `Isaac.jl`")
@testset "Isaac" begin

include("testproblems.jl")
include("test_nsoli.jl")
include("test_dlanczos.jl")
# include("test_saddle1.jl")
# include("test_minim.jl")


# if notCI
#     include("performance.jl")
# end

end
