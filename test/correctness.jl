
@testset "CTKSolvers" begin

@testset "dirder" begin
   println("testing `dirder`")
   for (f, df, init, randinit, name) in tests
      for i = 1:5
         x = randinit()
         w = rand(length(x)) - 0.5; w /= norm(w)
         err = norm( df(x) * w - CTKSolvers.dirder(x, w, f, f(x)) )
         @test err < 1e-6
      end
   end
end

@testset "dgmres" begin
   println("testing dgmres")
   for (f, df, init, randinit, name) in tests
      x = randinit()
      fx = f(x)
      for tol in (1e-1, 1e-2, 1e-3)
         maxiter = 100
         u, error, total_iters = dgmres(fx, f, x, tol, maxiter)
         err = norm(CTKSolvers.dirder(x, u, f, fx) - (-fx))
         @test err / norm(fx) < tol
         @test abs(error[end] - err) < 1e-7
         @test total_iters < maxiter
      end
   end
end

@testset "nsoli" begin
   println("testing nsoli")
   for (f, df, init, randinit, name) in tests
      for i = 1:5
         x0 = randinit()
         atol = 1e-5
         x, it_hist, ierr, x_hist = nsoli(x0, f, atol = atol, rtol = 0.0, maxit = 100)
         @test norm(f(x)) < atol
      end
   end
end

end
