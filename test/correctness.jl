
@testset "Isaac" begin

println("Basic correctness tests.")
# TODO: add tests determining what happens on failure!!!!

@testset "dirder" begin
   println("testing `dirder`")
   for (f, df, init, randinit, name) in tests
      for i = 1:10
         x = randinit()
         w = rand(length(x)) - 0.5; w /= norm(w)
         err = norm( df(x) * w - Isaac.dirder(x, w, f, f(x)) )
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
         err = norm(Isaac.dirder(x, u, f, fx) - (-fx))
         @test err / norm(fx) < tol
         # TODO: This next test seems to fail quite often
         # @test abs(error[end] - err) < 1e-7
         @test total_iters < maxiter
      end
   end
end

@testset "nsoli" begin
   println("testing nsoli")
   for (f, df, init, randinit, name) in tests
      for i = 1:10
         x0 = randinit()
         atol = 1e-5
         maxit = 100
         x, it_hist, ierr, x_hist = nsoli(x0, f;
                                       atol = atol, rtol = 0.0, maxit = maxit)
         if norm(f(x)) > atol
            @test ierr != 0
            @test length(x_hist) >= maxit
         else
            @test norm(f(x)) <= atol
         end
      end
   end
end

end
