using Isaac
using Base.Test

isCI = haskey(ENV, "CI")
notCI = !isCI

println("Running tests for `Isaac.jl`")
@testset "Isaac" begin

# include("testproblems.jl")
# include("correctness.jl")
# if notCI
#     include("performance.jl")
# end
# include("test_dlanczos.jl")
# include("test_saddle1.jl")
include("test_minim.jl")

end


# Two codes to look at for devising tests against NLopt
# https://github.com/JuliaDiffEq/DiffEqParamEstim.jl/blob/master/test/lorenz_true_test.jl
# https://github.com/JuliaDiffEq/DiffEqParamEstim.jl/blob/master/test/tests_on_odes/optim_test.jl#L16



# using SaddleSearch
# using SaddleSearch.TestSets
# using SaddleSearch.TestSets: hessprecond, precond
# using Isaac
# using Isaac: dlanczos, darnoldi, stabilise
#
#
# R = 5.1
# P = I
# V = LJVacancy2D(R=R, bc = :clamped)
# E, dE = objective(V)
# x0, v0 = ic_dimer(V, :near)
# P = precond(V, x0)
# Pprep = (P, x) -> precond(V, x)
#
# # test that the dlanczos and darnoldi return the same solutions
# # b = rand(length(x0))
# # xa, numfa, succa, isna = darnoldi( dE(x0), dE, x0, b, 1e-12, 10, Isaac.stabilise )
# # xl, _, numfl, succl, isnl = dlanczos( dE(x0), dE, x0, b, 1e-12, 10, Isaac.minim )
# # @show numfa, numfl, isna, isnl
# # @show norm(xa-xl, Inf)
#
# println("Testing nkminim, 2D LJ Vacancy, R = $(round(Int,R)), precon = $(!(P==I))")
# x, ndE = nkminim(E, dE, x0; verbose=0, P = P, precon_prep = Pprep)
# println("   num_dE = ", ndE)
# println("   |∇E| = ", norm(dE(x), Inf))
#
# println("Testing nsolistab, 2D LJ Vacancy, R = $(round(Int,R)), precon = $(!(P==I))")
# x, ndE = nsolistab(dE, x0; verbose=0, P = P, precon_prep = Pprep)
# println("   num_dE = ", ndE)
# println("   |∇E| = ", norm(dE(x), Inf))
