
using SaddleSearch
using SaddleSearch.TestSets
using SaddleSearch.TestSets: hessprecond, precond


println("Testing nsolistab for minimisation")
@testset "nsolidstab-minim" begin

#[1] p-Laplacian

plap(U; n=length(U)) = (n-1) * sum( (0.1 + diff(U).^2).^2 ) - sum(U) / (n-1)
plap1(U; n=length(U), dU = diff(U), dW = 4 * (0.1 + dU.^2) .* dU) =
                        (n-1) * ([0.0; dW] - [dW; 0.0]) - ones(U) / (n-1)
precond_plap(x::Vector) = precond_plap(length(x))
precond_plap(n::Number) = spdiagm( ( -ones(n-1), 2*ones(n), -ones(n-1) ),
                              (-1,0,1), n, n) * (n+1)
E_plap(X) = plap([0;X;0])
dE_plap(X) = plap1([0;X;0])[2:end-1]

for (n, prec) in [ (10, false), (10, true), (50, false), (50, true) ]
   x0 = zeros(n)
   P = prec ? precond_plap(n) : I
   println("Testing nsolistab, minim, p-Laplacian, n = $n, precon = $prec")
   x, ndE = nsolistab(dE_plap, x0, 0; krylovinit=:res, P = P)
   println("   num_dE = ", ndE)
   println("   |∇E| = ", norm(dE_plap(x), Inf))
   @test norm(dE_plap(x), Inf) < 1e-5
end


# [2] 2D LJ Vacancy Test

for (R, prec) in [(5.1, false), (5.1, true), (12.1, false), (12.1, true)]
   V = LJVacancy2D(R=R, bc = :clamped)
   E, dE = objective(V)
   x0, v0 = ic_dimer(V, :far)
   P, Pprep = prec ? (precond(V, x0), (P, x) -> precond(V, x)) : (I, (P,x)->P)

   println("Testing nsolistab, minim, 2D LJ Vacancy, R = $(round(Int,R)), precon = $(!(P==I))")
   x, ndE = nsolistab(dE, x0, 0; krylovinit=:res, P = P, precon_prep = Pprep)
   println("   num_dE = ", ndE)
   println("   |∇E| = ", norm(dE(x), Inf))
   @test norm(dE(x), Inf) < 1e-5
end


end  # @testset "nsolidstab-index1"
