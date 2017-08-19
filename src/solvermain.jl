using Parameters, ProgressMeter

export nsolimod, nkminim, nsolistab

"""
`nsolimod{T}(dE, x0::Vector{T}, saddleindex; kwargs...)`
   -> x, numdE

# A "modulated jacobian-free Newton-Krylov solver"

A Newton-Krylov solver for computing critical points of `E` with
prescribed `saddleindex`, e.g. `saddleindex = 0` returns only minima,
`saddleindex = 1` returns only index-1 saddles and so forth.
For computing *any* critical point use `nsoli` instead.

## Required Arguments

* `dE` : evaluate potential gradient
* `x0` : initial condition
* `saddleindex` : specifies which critical points are being sought, e.g.,
   `0` for minima, `1` for index-1 saddles; under some idealising assumptions,
   for `saddleindex=n`, the `nsolimod` dynamical system has as its stable equilibria all
   critical points with `n` negative and `d-n` positive hessian eigenvalues
   while all other critical points of `E` are unstable equilibria.

## Keyword Arguments

* `tol = 1e-5`
* `maxnumdE = 200`
* `maxstep = Inf`
* `hfd = 1e-7`
* `P = I, precon_prep = (P, x) -> P`
* `eigatol = 1e-1, eigrtol = 1e-1`
* `verbose = 1`
* `V0 = rand(T, (length(x0), saddleindex+1))`
* `E = nothing`
* `linesearch = nothing`
* `krylovinit = :resrot`

## Output

* `x` : solution
* `numdE` : number of gradient evaluations
"""
function nsolimod{T}(dE, x0::Vector{T}, saddleindex::Int;
                  tol = 1e-5,
                  maxnumdE = 200,
                  maxstep = Inf,
                  hfd = 1e-7,
                  P = I, precon_prep = (P, x) -> P,
                  eigatol = 1e-1, eigrtol = 1e-1,
                  verbose = 1,
                  V0 = rand(T, (length(x0), saddleindex)),
                  E = nothing,
                  linesearch = lsdefault,
                  krylovinit = :resrot,    # TODO: remove asap
                  update_α_old = true
               )
   debug = verbose > 2
   progressmeter = verbose == 1

   # initialise some more parameters; TODO: move these into NK?
   d = length(x0)
   eta = etamax = 0.5   # 0.9
   gamma = 0.9
   kmax = min(40, d)
   Carmijo = 1e-4
   α_old = α = 1.0
   σtransform = D -> indexp(D, saddleindex)

   # evaluate the initial residual
   x = copy(x0)
   v = copy(V0)
   f0 = dE(x)
   numdE = 1
   res = norm(f0, Inf)

   P = precon_prep(P, x)
   fnrm = nkdualnorm(P, f0)
   fnrmo = fnrm
   itc = 0

   while res > tol && numdE < maxnumdE
      rat = fnrm / fnrmo   # this is used for the Eisenstat-Walker thing
      itc += 1

      # compute the (modified) Newton direction
      if krylovinit == :res
         V0 = - (P \ f0)
      elseif krylovinit == :rand
         V0 = P \ rand(d)
      elseif krylovinit == :rot
         V0 = v
      elseif krylovinit == :resrot
         V0 = [ P\f0  v ]
      else
         error("unknown parameter `krylovinit = $(krylovinit)`")
      end

      p, G, inner_numdE, success, isnewton =
            dlanczos(f0, dE, x, -f0, eta * nkdualnorm(P, f0), kmax, σtransform;
                         P = P, V0 = V0, debug = debug,
                         hfd = hfd, eigatol = eigatol, eigrtol = eigrtol)

      numdE += inner_numdE
      if debug; @show isnewton; end

      # ~~~~~~~~~~~~~~~~~~ LINESEARCH ~~~~~~~~~~~~~~~~~~~~~~
      # output of linesearch will be: α, xt (new x), ft (new dE), nft (norm)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if isnewton   # if we are in the Newton regime, then we always try to reduce the residual
         iarm = 0
         α = αt = 1.0
         xt = x + αt * p
         ft = dE(xt)
         numdE += 1
         nf0 = nkdualnorm(P, f0)
         nft = nkdualnorm(P, ft)

         αm = 1.0  # these have no meaning; just allocations
         nfm = 0.0

         while nft > (1 - Carmijo * αt) * nf0
            if iarm == 0
               α *= 0.5
            else
               α = parab3p( αt, αm, nf0^2, nft^2, nfm^2; sigma0 = 0.1, sigma1 = 0.5 )
            end
            αt, αm, nfm = α, αt, nft
            if αt < 1e-8   # TODO: make this a parameter
               error(" Armijo failure, step-size too small")
            end

            xt = x + αt * p
            ft = dE(xt)
            numdE += 1
            nft = nkdualnorm(P, ft)
            iarm += 1
         end
         α_old = αt
      else
         # if we are here, then p is not a newton direction
         # (i.e. an e-val was flipped)
         αt, xt, ft, nft, numdE_plus = linesearch(x, p, α_old, E, dE, f0, P, maxstep, G)
         numdE += numdE_plus
         # TODO: we get better results if we DO NOT update α_old here.
         #       WHY?????
         if update_α_old
            α_old = αt
         end
      end
      if debug; @show αt; end
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      # update current configuration and preconditioner
      x, f0, fnrm = xt, ft, nft
      P = precon_prep(P, x)
      res = norm(f0, Inf)
      fnrm = nkdualnorm(P, f0)     # should be the dual norm!
      rat = fnrm/fnrmo

      # if debug; @show λ, res; end

      if res <= tol
         return x, numdE
      end

      # ---------------------------------------------------------
      # Adjust eta as per Eisenstat-Walker.   # TODO: make this a flag!
      # TODO: check also that we are in the index-1 regime (what do we do if
      # not? probably reset eta?)
      # S. C. Eisenstat, H. F. Walker,
      # Choosing the forcing terms in an inexact Newton method,
      # SIAM J. Sci. Comput. 17 (1996) 16–32.
      etaold = eta
      etanew = gamma * rat^2
      if gamma * etaold^2 > 0.1
         etanew = max(etanew, gamma * etaold^2)
      end
      eta = min(etanew, etamax)
      eta = max(eta, 0.5 * tol / fnrm)
      # ---------------------------------------------------------
   end

   if verbose >= 1
      warn("NK did not converge within the maximum number of dE evaluations")
   end
   return x, numdE
end



function nkminim(E, dE, x0;
                tol = 1e-5,
                maxnumdE = 200,
                maxstep = Inf,
                hfd = 1e-7,
                P = I, precon_prep = (P, x) -> P,
                verbose = 1,
                update_α_old = true,
                linesearch = lsarmijo )

   # specialised settings
   eigatol = Inf
   eigrtol = Inf
   krylovinit = :res

   # call the generic solver
   return nsolimod(dE, x0, 0;
            tol = tol, maxnumdE = maxnumdE, maxstep = maxstep, hfd = hfd,
            P = P, precon_prep = precon_prep, eigatol = eigatol, eigrtol = eigrtol,
            verbose = verbose, linesearch = linesearch, krylovinit = krylovinit,
            E = E, update_α_old = update_α_old)
end






"""
`nsolistab{T}(f, x0::Vector{T}; kwargs...) -> x, numdE`

# A "modulated jacobian-free Newton-Krylov solver"

A Newton-Krylov solver for solving stable equilibria of ẋ = -f(x), i.e.,
solutions to f(x) = 0 satisfying Reλ > 0 for all λ \in σ(∂f(x)). The choice
of -f as opposed to f is motivated by the gradient flow case.

## Required Arguments

* `f` : evaluate nonlinear system
* `x0` : initial condition

## Keyword Arguments

* `tol = 1e-5`
* `maxnumdE = 200`
* `maxstep = Inf`
* `hfd = 1e-7`
* `P = I, precon_prep = (P, x) -> P`
* `verbose = 1`
* `linesearch = lsdefault`

## Output

* `x` : solution
* `numdE` : number of gradient evaluations
"""
function nsolistab{T}(f, x0::Vector{T};
                  tol = 1e-5,
                  maxnumdE = 200,
                  maxstep = Inf,
                  hfd = 1e-7,
                  P = I, precon_prep = (P, x) -> P,
                  verbose = 1,
                  linesearch = lsdefault )
   debug = verbose > 2
   progressmeter = verbose == 1

   # initialise some more parameters
   d = length(x0)
   eta = etamax = 0.5   # 0.9 ???
   gamma = 0.9
   kmax = min(40, d)
   Carmijo = 1e-4
   α_old = α = 1.0

   # evaluate the initial residual
   x = copy(x0)
   f0 = f(x)
   numdE = 1
   res = norm(f0, Inf)

   P = precon_prep(P, x)
   fnrm = norm(f0)
   fnrmo = fnrm
   itc = 0

   while res > tol && numdE < maxnumdE
      rat = fnrm / fnrmo   # this is used for the Eisenstat-Walker thing
      itc += 1

      p, inner_numdE, success, isnewton =
            darnoldi(f0, f, x, -f0, eta * norm(f0), kmax, stabilise;
                         P = P, debug = debug, hfd = hfd)

      numdE += inner_numdE
      if debug; @show isnewton; end

      # ~~~~~~~~~~~~~~~~~~ LINESEARCH ~~~~~~~~~~~~~~~~~~~~~~
      # output of linesearch will be: α, xt (new x), ft (new dE), nft (norm)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if isnewton   # if we are in the Newton regime, then we always try to reduce the residual
         iarm = 0
         α = αt = 1.0
         xt = x + αt * p
         ft = f(xt)
         numdE += 1
         nf0 = norm(f0)
         nft = norm(ft)

         αm = 1.0  # these have no meaning; just allocations
         nfm = 0.0

         while nft > (1 - Carmijo * αt) * nf0
            if iarm == 0
               α *= 0.5
            else
               α = parab3p( αt, αm, nf0^2, nft^2, nfm^2; sigma0 = 0.1, sigma1 = 0.5 )
            end
            αt, αm, nfm = α, αt, nft
            if αt < 1e-8   # TODO: make this a parameter
               error(" Armijo failure, step-size too small")
            end

            xt = x + αt * p
            ft = f(xt)
            numdE += 1
            nft = norm(ft)
            iarm += 1
         end
         α_old = αt
      else
         # if we are here, then p is not a newton direction (i.e. an e-val was flipped)
         αt, xt, ft, nft, numdE_plus = linesearch(x, p, α_old, nothing, f, f0, I, maxstep, I)
         numdE += numdE_plus
         α_old = αt
      end
      if debug; @show αt; end
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      # update current configuration and preconditioner
      x, f0, fnrm = xt, ft, nft
      P = precon_prep(P, x)
      res = norm(f0, Inf)
      fnrm = norm(f0)     # should be the dual norm!
      rat = fnrm/fnrmo

      if debug; @show res; end

      if res <= tol
         return x, numdE
      end

      # ---------------------------------------------------------
      # Adjust eta as per Eisenstat-Walker.   # TODO: make this a flag!
      # TODO: check also that we are in the newton regime (what do we do if
      # not? probably reset eta?)
      # S. C. Eisenstat, H. F. Walker,
      # Choosing the forcing terms in an inexact Newton method,
      # SIAM J. Sci. Comput. 17 (1996) 16–32.
      etaold = eta
      etanew = gamma * rat^2
      if gamma * etaold^2 > 0.1
         etanew = max(etanew, gamma * etaold^2)
      end
      eta = min(etanew, etamax)
      eta = max(eta, 0.5 * tol / fnrm)
      # ---------------------------------------------------------
   end

   if verbose >= 1
      warn("NK did not converge within the maximum number of dE evaluations")
   end
   return x, numdE
end
