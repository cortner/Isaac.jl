
# TODO: replace with with LineSearches or alternatively move into LineSearches

"""
`function  parab3p(lambdac, lambdam, ff0, ffc, ffm) -> lambdap`

Apply three-point safeguarded parabolic model for a line search.

### Input Parameters

* `lambdac` : current steplength
* `lambdam` : previous steplength
* `ff0` : value of |F(x_c)|^2
* `ffc` : value of |F(x_c + λc d)|^2
* `ffm` : value of |F(x_c + λm d)|^2

### Output

* lambdap : new value of lambda given parabolic model

### Keyword arguments (internal parameters)

* `sigma0 = .1`, `sigma1 = .5` : safeguarding bounds for the linesearch
"""
function parab3p(lambdac, lambdam, ff0, ffc, ffm; sigma0 = 0.1, sigma1 = 0.5)
   # compute coefficients of interpolation polynomial
   # p(lambda) = ff0 + (c1 lambda + c2 lambda^2)/d1
   # d1 = (lambdac - lambdam)*lambdac*lambdam < 0
   #      so if c2 > 0 we have negative curvature and default to
   #      lambdap = sigam1 * lambda
   c2 = lambdam * (ffc - ff0) - lambdac * (ffm - ff0)
   if c2 >= 0
      return sigma1 * lambdac
   end
   c1 = lambdac^2 * (ffm-ff0) - lambdam^2 * (ffc - ff0)
   lambdap = -c1 * 0.5 / c2
   if lambdap < sigma0 * lambdac
      return sigma0 * lambdac
   elseif lambdap > sigma1 * lambdac
      return sigma1 * lambdac
   end
   return lambdap
end



function parab2p(f0, df0, f1, alpha)
   return - (df0 * alpha^2) / ( 2.0 * (f1 - f0 - df0*alpha) )
end

"""
at the moment, this is just quadratic back-tracking

TODO:
 * probably switch to LineSearches.jl
"""
function lsarmijo(x, p, α_old, E, dE, dE0, P, maxstep, K = nothing;
                  Ca = 0.2)

   maxα = maxstep / norm(p, Inf)
   # take one potentially forward step
   slope = dot(dE0, p)
   E0 = E(x)
   #  α1 = 0.5 + 0.33*α_old : works really  well for P = I
   #  α1 = 0.66*α_old : works really well with a good P
   #  α1 = 1.0 : decent overall
   if P == I
      α1 = 0.5 + 0.33 * α_old
   else
      α1 = 0.66 * α_old
   end
   α1 = min(α1, maxα)
   x1 = x + α1 * p
   E1 = E(x1)
   numdE = 1

   if E1 <= E0 + Ca * α1 * slope
      # Armijo is satisfied, we will now try to improve a bit but won't
      # try too hard and just revert to x1 if needed
      if E0 + α1 * slope < E1   # above the slope > parab2p is meaningful
         αt = parab2p(E0, slope, E1, α1)
      else # below the slope > parab2p would give nonsense
         αt = 1.5 * max(α1, α_old)
      end
      # make sure the step is not too small
      αt = max(α1/4, αt)
      # and not too large
      αt = min(α1*2, maxα, αt)
      # update configuration
      xt = x + αt * p
      Et = E(xt)
      numdE += 1
      # check whether to return xt or x1
      if Et > E1
         Et, αt, xt = E1, α1, x1
      end
   else
      # CASE 2: the starting guess does not satisfy Armijo
      # > start a back-tracking loop
      xt, αt, Et = x1, α1, E1
      while Et > E0 + Ca * αt * slope
         αt_old, αt = αt, parab2p(E0, slope, Et, αt)
         αt = max(αt_old/4, αt)
         xt = x + αt * p
         Et = E(xt)
         numdE += 1
      end
   end

   # we also need to return the gradient at the final l.s. point
   dEt = dE(xt)
   nft = nkdualnorm(P, dEt)
   return αt, xt, dE(xt), nft, numdE
end


"""
A very crude step-length selection for the pre-asymptotic regime
"""
function lsforce(x, p, α_old, E, dE, f0, P, maxstep, K;
            Ca = 0.2, Cw = 0.9, minα = 1e-6, αinit = :one)
   # in this case, we do something very crude:
   #   take same step as before, then do one line-search step
   #   and pick the better of the two.
   φ = f_ -> dot(K \ f_, P, p)  # dot(f_, p)   # TODO: not sure which one will be better
   g0 = φ(f0)
   tol = Cw * abs(g0)
   maxα = maxstep / norm(p, Inf)
   numdE = 0
   if αinit == :one
      α1 = 1.0
   elseif αinit == :heuristic
      if P == I
         α1 = 0.5 + 0.33 * α_old
      else
         α1 = 0.66 * α_old
      end
   else
      error("unknown `αinit`")
   end
   α1 = min(α1, maxα)
   x1 = x + α1 * p
   f1 = dE(x1)
   g1 = φ(f1)
   numdE += 1

   acceptstep = gg -> 2 * g0 < gg < 0.2 * abs(g0)

   if acceptstep(g1)
      # put a straight line through (0, g0), (α1, g1) and intersect with 0
      if abs(g1-g0)/α1 < 1e-3
         αt = 0.5 * α1
      else
         αt = - g0 * α1 / (g1-g0)   # g0 + α * (g1-g0)/α1 = 0
      end
      # make sure the step is not too small
      αt = max(α1/4, αt)
      # and not too large
      αt = min(α1*2, maxα, αt)

      # update configuration
      xt = x + αt * p
      ft = dE(xt)
      gt = φ(ft)
      numdE += 1
      # check whether to return xt or x1
      if !acceptstep(gt) || abs(gt) > abs(g1)
         gt, αt, xt = g1, α1, x1
      end
   else
      αt, gt, ft, xt = α1, g1, f1, x1
      # reduce step-length until the new gradient is not too bad
      while (2 * g0 > gt) || (gt > 0.2 * abs(g0))
         αt *= 0.5
         xt = x + αt * p
         ft = dE(xt)
         numdE += 1
         gt = φ(ft)
      end
   end

   return αt, xt, ft, dualnorm(P, ft), numdE
end
