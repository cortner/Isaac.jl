using Parameters


"""
Newton-Krylov based saddle search method
"""
@with_kw type NK
   tol::Float64 = 1e-5
   maxnumdE::Int = 1000
   len::Float64 = 1e-7
   precon = I
   precon_prep! = (P, x) -> P
   verbose::Int = 1
   krylovinit::Symbol = :resrot  # allow res, rand, rot, resrot
   maxstep::Float64 = Inf
   eigatol::Float64 = 1e-1
   eigrtol::Float64 = 1e-1
end


function run!{T}(method::NK, E, dE, x0::Vector{T},
                  v0::Vector{T} = rand(T, length(x0)) )
   # get parameters
   @unpack tol, maxnumdE, len, verbose, krylovinit, maxstep,
      eigatol, eigrtol = method
   precon = x -> method.precon_prep!(method.precon, x)
   debug = verbose > 2

   # initialise some more parameters; TODO: move these into NK?
   d = length(x0)
   eta = etamax = 0.5   # 0.9
   gamma = 0.9
   kmax = min(40, d)
   Carmijo = 1e-4
   α_old = α = 1.0

   # evaluate the initial residual
   x = copy(x0)
   v = copy(v0)
   f0 = dE(x)
   numdE = 1
   res = norm(f0, Inf)

   P = precon(x)
   fnrm = nkdualnorm(P, f0)
   fnrmo = 1.0
   itc = 0

   while res > tol && numdE < maxnumdE
      rat = fnrm / fnrmo   # TODO: is this unused?
      fnrmo = fnrm         # TODO: probably move this to where fnrm is updated!
      itc += 1

      if debug; @show dot(f0, v); end
      # compute the (modified) Newton direction
      if krylovinit == :res
         V0 = - P \ f0
      elseif krylovinit == :rand
         V0 = P \ rand(d)
      elseif krylovinit == :rot
         V0 = v
      elseif krylovinit == :resrot
         V0 = [ P \ f0 v ]
      else
         error("unknown parameter `krylovinit = $(krylovinit)`")
      end
      p, λ, v, inner_numdE, isnewton =
            blocklanczos(f0, dE, x, eta * norm(f0), kmax;
                         P = P, b = - f0, V0 = V0, debug = (verbose >= 3),
                         h = len, eigatol = eigatol, eigrtol = eigrtol)
      numdE += inner_numdE
      if debug; @show isnewton; end

      # ~~~~~~~~~~~~~~~~~~ LINESEARCH ~~~~~~~~~~~~~~~~~~~~~~
      if isnewton
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
         # if we are here, then p is not a newton direction (i.e. an e-val was
         # flipped) in this case, we do something very crude:
         #   take same step as before, then do one line-search step
         #   and pick the better of the two.
         αt = 0.66 * α_old  # probably can do better by re-using information from dcg_...
         αt = min(αt, maxstep / norm(p, Inf))
         xt = x + αt * p
         ft = dE(xt)
         nft = nkdualnorm(P, ft)
         numdE += 1
         # find a root: g(t) = (1-t) f0⋅p + t ft ⋅ p = 0 => t = f0⋅p / (f0-ft)⋅p
         #    if (f0-ft)⋅p is very small, then simply use xt as the next step.
         #    if it is large enough, then take just one iteration to get a root
         if abs(dot(f0 - ft, P, p)) > 1e-4   #  TODO: make this a parameter
            t = dot(f0, P, p) / dot(f0 - ft, P, p)
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
      end
      if verbose > 3; @show αt; end
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      # update current configuration and preconditioner
      x, f0, fnrm = xt, ft, nft
      P = precon(x)
      res = norm(f0, Inf)
      fnrm = nkdualnorm(P, f0)     # should be the dual norm!
      rat = fnrm/fnrmo

      if verbose > 3; @show λ, res; end

      if res <= tol
         return x, numdE
      end

      # ---------------------------------------------------------
      # Adjust eta as per Eisenstat-Walker.   # TODO: make this a flag!
      # TODO: check also the we are in the index-1 regime (what do we do if
      # not? probably reset eta?)
      etaold = eta
      etanew = gamma * rat^2
      if gamma * etaold^2 > 0.1
         etanew = max(etanew, gamma * etaold^2)
      end
      eta = min(etanew, etamax)
      eta = max(eta, 0.5 * tol / fnrm)
      # ---------------------------------------------------------

   end

   if verbose > 1
      warn("NK did not converge within the maximum number of dE evaluations")
   end
   return x, numdE
end
