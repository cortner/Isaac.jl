
atol = 1e-5
rtol = 0.0
maxit = 100

println("=====================")
println("  Performance Tests ")
println("=====================")
for (f, df, init, _, name) in tests
   println("-------------------")
   println("  $name Test")
   println("-------------------")
   x, it_hist, ierr, x_hist = nsoli(init(), f,
                                 atol = atol, rtol = rtol, maxit = maxit)
   # improve the solution a bit
   y = copy(x)
   for i = 1:3
      y -= df(y) \ f(y)
   end
   println("     Error Message: ", ierr)
   println("      # iterations: ", length(x_hist) )
   println("          Residual: ", norm(f(x), Inf))
   println("     Residual xref: ", norm(f(y), Inf))
   println("  Error |x - xref|: ", norm(x - y, Inf))
end
