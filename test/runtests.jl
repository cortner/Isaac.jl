using CTKSolvers
using ForwardDiff
using Base.Test


f(x) = [x[1] + x[2] + x[1]^2, x[2] + x[1]*x[2]]
df(x) = ForwardDiff.jacobian(f, x)

@testset "CTKSolvers" begin

@testset "dirder" begin
   println("testing `dirder`")
   for n = 1:10
      x = rand(2)
      w = rand(2); w /= norm(w)
      err = norm( df(x) * w - CTKSolvers.dirder(x, w, f, f(x)) )
      @test err < 1e-6
   end
end

@testset "dgmres" begin
   println("testing dgmres")
   x = rand(2)
   fx = f(x)
   for tol in (1e-1, 1e-2, 1e-3)
      maxiter = 100
      u, error, total_iters = dgmres(fx, f, x, tol, maxiter)
      err = norm(CTKSolvers.dirder(x, u, f, fx) - (-fx))
      @test err < tol
      @test abs(error[end] - err) < 1e-7
      @test total_iters < maxiter
   end
end

@testset "nsoli" begin
   println("testing nsoli")
   for n = 1:10
      x0 = rand(2)
      atol = 1e-5
      x, it_hist, ierr, x_hist = nsoli(x0, f, atol = atol, rtol = 0.0, maxit = 100)
      @test norm(f(x)) < atol
   end
end

end
