
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
   #      lambdap = sigma1 * lambda
   c2 = lambdam * (ffc - ff0) - lambdac * (ffm - ff0)
   if c2 >= 0
      return sigma1 * lambdac
   end
   c1 = lambdac^2 * (ffm - ff0) - lambdam^2 * (ffc - ff0)
   lambdap = -c1 * 0.5 / c2
   if lambdap < sigma0 * lambdac
      return sigma0 * lambdac
   elseif lambdap > sigma1 * lambdac
      return sigma1 * lambdac
   end
   return lambdap
end
