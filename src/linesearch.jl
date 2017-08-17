
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

"""
A very crude step-length selection for the pre-asymptotic regime.

TODO: remove P from the arguments!
"""
function onestep(x, p, α_old, E, dE, f0, P, maxstep)
   # in this case, we do something very crude:
   #   take same step as before, then do one line-search step
   #   and pick the better of the two.
   numdE = 0
   αt = 0.66 * α_old
   αt = min(αt, maxstep / norm(p, Inf))
   xt = x + αt * p
   ft = dE(xt)
   numdE += 1
   nft = nkdualnorm(P, ft)
   # find a root: g(t) = (1-t) f0⋅p + t ft ⋅ p = 0 => t = f0⋅p / (f0-ft)⋅p
   #    if (f0-ft)⋅p is very small, then simply use xt as the next step.
   #    if it is large enough, then take just one iteration to get a root
   if abs(dot(f0 - ft, P, p)) > 1e-4   #  TODO: make this a parameter
                                       #  TODO: there should be no P here!
      t = dot(f0, P, p) / dot(f0 - ft, P, p)     # . . . nor here!
      t = max(t, 0.1)    # don't make too small a step
      t = min(t, 4 * t)  # don't make too large a step
      αt, αm, nfm, fm = (t*αt), αt, nft, ft
      αt = min(αt, maxstep / norm(p, Inf))
      xt = x + αt * p
      ft = dE(xt)
      numdE += 1
      nft = nkdualnorm(P, ft)
      # if the line-search step is worse than the initial trial, then
      # we revert
      if abs(dot(ft, P, p)) > abs(dot(fm, P, p))
         αt, xt, nft, ft = αm, x + αm * p, nfm, fm
      end
   end
   return αt, xt, ft, nft, numdE
end


function parab2p(f0, df0, f1, alpha)
   return - (df0 * alpha^2) / ( 2.0 * (f1 - f0 - df0*alpha) )
end

"""
at the moment, this is just quadratic back-tracking

TODO:
 * switch to LineSearches.jl
 *
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

@noinline function _dotK_(p, f, P, K)
   
   fK = K \ f
   return dot(p, P, fK)
end


"""
A very crude step-length selection for the pre-asymptotic regime.

TODO: remove P from the arguments!
"""
function lswolfe(x, p, α_old, E, dE, f0, P, maxstep, K; Cw = 0.9, minα = 1e-6)
   # in this case, we do something very crude:
   #   take same step as before, then do one line-search step
   #   and pick the better of the two.
   φ = f_ -> _dotK_(p, f_, P, K)
   g0 = φ(f0)
   tol = Cw * abs(g0)
   maxα = maxstep / norm(p, Inf)
   numdE = 0
   if P == I
      α1 = 0.5 + 0.33 * α_old
   else
      α1 = 0.66 * α_old
   end
   α1 = min(α1, maxα)
   x1 = x + α1 * p
   f1 = dE(x1)
   g1 = φ(f1)
   numdE += 1

   # check whether initial trial satisfies the Wolfe condition
   if abs(g1) <= tol
      # if yes, then we just take one extra step to try and find a better
      # guess, but if we cannot then we continue with αt

      # find a root: g(t) = (1-t) f0⋅p + t ft ⋅ p = 0 => t = f0⋅p / (f0-ft)⋅p
      #    if (f0-ft)⋅p is very small, then simply use xt as the next step.
      #    if it is large enough, then take just one iteration to get a root
      t = g0 / (g0-g1)
      t = max(t, 0.25)    # don't make too small a step
      t = min(t, 4.0)  # don't make too large a step
      αt = min(t * α1, maxα)
      xt = x + αt * p
      ft = dE(xt)
      gt = φ(ft)
      numdE += 1
      # if the line-search step is worse than the initial trial, then
      # we revert
      if abs(gt) > abs(g1)
         αt, xt = α1, x1
      end
   else
      # CASE 2: the initial guess does not satisfy Wolfe
      #         >>> need to work a bit now, hopefully it is worth it
      # let's try basic secant for now
      A = [0.0, α1]
      X = [x, x1]
      F = [f0, f1]
      G = [g0, g1]
      @show α1, g0, g1
      for n = 1:5  # dont try for too long
         αt = A[end] - G[end] * (A[end] - A[end-1])/(G[end]-G[end-1])
         αt = min(αt, maxα)
         αt = max(αt, minα)
         xt = x + αt * p
         ft = dE(xt)
         numdE += 1
         gt = φ(ft)
         @show αt, gt
         push!(A, αt)
         push!(X, xt)
         push!(F, ft)
         push!(G, gt)
         if (abs(gt) <= tol) || (αt == maxα) || (αt == minα)
            break
         end
      end
      # now pick the best step we have found
      A = A[2:end]
      X = X[2:end]
      F = F[2:end]
      G = G[2:end]
      _, imin = findmin(abs.(G))
      αt, xt, ft = A[imin], X[imin], F[imin]
   end

   return αt, xt, ft, dualnorm(P, ft), numdE
end
