
"""
function nsoli(x,f; kwargs...) -> (sol, it_hist, ierr, x_hist)

Newton-Krylov solver, globally convergent solver for f(x) = 0, using
Inexact-Newton-Armijo iteration, Eisenstat-Walker forcing term
and Parabolic line search via three point interpolation.

This code has been transcribed with permission from
http://www4.ncsu.edu/~ctk/newton/SOLVERS/nsoli.m
C. T. Kelley, April 27, 2001

### Required Input Parameters

* x: initial iterate
* f: function for which we are seeking a root

### Output Parameters

* sol = solution
* it_hist(maxit,3) = l2 norms of nonlinear residuals
           for the iteration, number of function evaluations,
           and number of steplength reductions
* ierr: 0 upon successful termination; 1 if after maxit iterations
           the termination criterion is not satsified; 2 failure in the line search. The iteration
           is terminated if too many steplength reductions
           are taken.
* x_hist: matrix of the entire interation history.
           The columns are the nonlinear iterates. This
           is useful for making movies, for example, but
           can consume way too much storage. This is an
           OPTIONAL argument. Storage is only allocated
           if x_hist is in the output argument list.

### Keyword Parameters

* atol: absolute error tolerance for the nonlinear iteration
* rtol: relative error tolerance for the nonlinear iteration
* lmeth: choice of linear iterative method
        1 (GMRES, no restarts), 2 GMRES(m), 3 (BICGSTAB), 4 (TFQMR);
     default = 1
* restart_limit: max number of restarts for GMRES if lmeth = 2, default = 20
* maxit: maxmium number of nonlinear iterations, default = 40
* maxitl: maximum number of inner iterations before restart in GMRES(m), m = maxitl, default = 40. For iterative methods other than GMRES(m) maxitl is the upper bound on linear iterations.
* |etamax| : Maximum error tolerance for residual in inner
     iteration. The inner iteration terminates
     when the relative linear residual is
     smaller than eta*| F(x_c) |. eta is determined
     by the modified Eisenstat-Walker formula if etamax > 0.
     If etamax < 0, then eta = |etamax| for the entire
     iteration. default: etamax = .9
* debug = turns on/off iteration statistics display as
               the iteration progresses
* alpha = 1.d-4, parameter to measure sufficient decrease (Armijo)
* sigma0 = .1, sigma1 = .5, safeguarding bounds for the linesearch
* maxarm = 20, maximum number of steplength reductions before
                    failure is reported
* debug = 0:  Set the debug parameter; 1 turns display on, otherwise off.
"""
function nsoli{T}(x::Vector{T}, f;
   atol = 1e-5, rtol = 1e-5, dkrylov = dgmres,
   maxit = 40, lmaxit = 40, etamax = 0.9, reorth = 1,
   alpha = 1e-4, sigma0 = 0.1, sigma1 = 0.5, maxarm = 20, gamma = 0.9,
   debug = 0 )

   # there was a paremter restart_limit = 20,  >>>> move into defn of dkrylov

   # Initialize some variables
   ierr = 0
   it_histx = zeros(maxit,3)
   x_hist = Vector{T}[x]

   # Initialize parameters for the iterative methods, Check for optional inputs.
   it_histx = zeros(maxit+1, 3)
   # gmparms = [abs(etamax), lmaxit, parms(5), 1]; TODO: remove this
   eta = abs(etamax)
   n = length(x)
   itc = 0          # iteration counter

   # Evaluate f at the initial iterate, and compute the stop tolerance.
   f0 = f(x)                             # TODO: count f evaluations instead of iterations?
   fnrm = norm(f0)                       # TODO: switch to Inf-norm???
   it_histx[itc+1, :] = [fnrm, 0, 0]

   fnrmo = 1.0      # old fnrm
   stop_tol = atol + rtol * fnrm

   # TODO: replace this with a nice type collecting iteration information
   outstat = zeros(1, 5)
   outstat[itc+1, :] = [itc, fnrm, 0, 0, 0]

   # main iteration loop
   while fnrm > stop_tol && itc < maxit
      # Keep track of the ratio (rat = fnrm/frnmo) of successive residual norms
      # and the iteration counter (itc).
      rat = fnrm / fnrmo
      fnrmo = fnrm
      itc = itc + 1

      # compute the Newton direction # TODO: is zero the best initialisation?
      step, errstep, inner_it_count, inner_f_evals =
            dkrylov(f0, f, x, eta, lmaxit, reorth)

      # ~~~~~~~~~~~~~~ BEGINNING OF LINESEARCH ~~~~~~~~~~~~~~~~
      xold = x
      lambda = 1
      lamm = 1
      lamc = lambda
      iarm = 0
      xt = x + lambda * step
      ft = f(xt)
      nft = norm(ft)
      nf0 = norm(f0)
      ff0 = nf0 * nf0
      ffc = nft * nft
      ffm = nft * nft
      while nft >= (1 - alpha * lambda) * nf0    # |ft| >= (1 - α λ) |f0|
         # Apply the three point parabolic model.
         if iarm == 0
            lambda = sigma1 * lambda
         else
            lambda = parab3p( lamc, lamm, ff0, ffc, ffm;
                              sigma0 = sigma0, sigma1 = sigma1 )
         end

         # Update x; keep the books on lambda.
         xt = x + lambda * step
         lamm = lamc
         lamc = lambda

         # Keep the books on the function norms.
         ft = f(xt)
         nft = norm(ft)
         ffm = ffc
         ffc = nft*nft
         iarm = iarm+1
         if iarm > maxarm
            warn(" Armijo failure, too many reductions ")
            ierr = 2
            display(outstat)
            it_hist = it_histx[1:itc+1, :]
            push!(x_hist, x)
            sol = xold
            return sol, it_hist, ierr, x_hist
         end
      end
      x = xt
      f0 = ft
      # ~~~~~~~~~~~~~~~ END OF LINESEARCH ~~~~~~~~~~~~~~~~

      push!(x_hist, x)
      fnrm = norm(f0)
      it_histx[itc+1, 1] = fnrm     # TODO: there is a bug here when we reach too many iterations

      # How many function evaluations did this iteration require?
      it_histx[itc+1, 2] = it_histx[itc, 2] + inner_f_evals + iarm + 1
      if itc == 1
         it_histx[itc+1, 2] = it_histx[itc+1, 2] + 1
      end
      it_histx[itc+1, 3] = iarm

      rat = fnrm/fnrmo
      # Adjust eta as per Eisenstat-Walker.   # TODO: make this a flag!
      if etamax > 0
         etaold = eta
         etanew = gamma * rat^2
         if gamma * etaold^2 > 0.1
            etanew = max(etanew, gamma * etaold^2)
         end
         eta = min(etanew, etamax)
         eta = max(eta, 0.5 * stop_tol / fnrm)
      end

      outstat = [outstat; [itc fnrm inner_it_count rat iarm]]
      #  outstat[itc+1, :] = [itc fnrm inner_it_count rat iarm]
   end
   sol = x
   it_hist = it_histx[1:itc+1, :]
   if debug == 1
      display(outstat)
      it_hist = it_histx[1:itc+1, :]
   end

   # on failure, set the error flag
   if fnrm > stop_tol
      ierr = 1
   end

   return sol, it_hist, ierr, x_hist
end
