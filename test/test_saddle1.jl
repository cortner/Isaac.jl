
using SaddleSearch
using SaddleSearch.TestSets
using SaddleSearch: numE, numdE, res_trans, res_rot
using SaddleSearch.TestSets: hessprecond, precond
using Isaac
using SaddleSearch: ODE12r, odesolve, IterationLog

println("Testing nsolimod for index-1 saddles")
@testset "nsolidmod-index1" begin

#[1] Muller Test
xe = []
for init in (:near, :far)
   println("Testing nsolimod, Muller, init = $(init)")
   V = MullerPotential()
   x0, v0 = ic_dimer(V, init)
   E, dE = objective(V)
   x, ndE = nsolimod(dE, x0, 1; V0=v0)
   println("   num_dE = ", ndE)
   if init == :near
      xe = copy(x)
   else
      println("  |x - xe| = ", norm(x - xe, Inf))
      @test norm(x - xe, Inf) < 1e-8
   end
end

# [2] Vacancy Test
V = LJVacancy2D(bc = :clamped)
E, dE = objective(V)
x0, v0 = ic_dimer(V, :near)
P = precond(V, x0)
Iprep = (P, x) -> P
Pprep = (P, x) -> precond(V, x)
xe = []

for (init, precond, precon_prep) in
      [(:near, I, Iprep), (:near, P, Pprep), (:far, I, Iprep), (:far, P, Pprep)]

   x0, v0 = ic_dimer(V, init)
   pre = !(precond == I)

   println("Testing nsolimod, 2D LJ Vacancy, init = $init, precond = $pre")
   x, ndE = nsolimod(dE, x0, 1; V0=v0, P = precond, precon_prep=precon_prep)
   println("   num_dE = ", ndE)
   @test norm(dE(x), Inf) < 1e-5
   println("   |dE| = ", norm(dE(x), Inf))
   if xe == []
      # to make sure all 4 methods converge to the same saddle,
      # store the first one ...
      xe = x
      H2 = TestSets.hessian(V, x)
      σ = sort(eigvals(Symmetric(H2)))
      @test σ[1] * σ[2] < 0
      println("   min-eigs = ", σ[1:2])
   else
      # ... and then compare the others
      @test norm(x - xe, Inf) < 1e-6
      println("   |x - xe| = ", norm(x - xe, Inf))
   end
end



end  # @testset "nsolidstab-index1"
