

export dgmres, dlanczos


"""
`dirder(x,w,f,f0)`

Finite difference directional derivative, approximate f'(x) w

### Inputs
* x, w = point and direction
* f = function
* f0 = f(x), in nonlinear iterations f(x) has usually been computed
                before the call to dirder

### Output
* `z`: approximation to `f'(x) w`

based on a code by C. T. Kelley, April 1, 2003
"""
function dirder(x, w, f, f0, h)
   # initialise difference increment
   epsnew = h
   n = length(x)
   # scale the step
   if norm(w) == 0
      return zeros(n)
   end
   # Now scale the difference increment.
   xs = dot(x, w) / norm(w)
   if xs != 0.0
      epsnew *= max(abs(xs), 1.0) * sign(xs)
   end
   epsnew /= norm(w)
   f1 = f(x + epsnew * w)
   return (f1 - f0) / epsnew
end


"""
`dirderinf(f, f0, xc, v, h)`:
directional derivative operator to mimic the ∂f(xc) * v action
via (f(x+h'*v) - f(x)) / h' where h' is rescaled against ||v||_∞
"""
function dirderinf(f, f0, xc, v, h)
   h = h / norm(v, Inf)
   return (f(xc + h * v) - f0) / h
end



"""
`function givapp(c, s, vin, k) -> vrot`

Apply a sequence of k Givens rotations, used within gmres codes.
"""
function givapp(c, s, vin, k)
   vrot = copy(vin)                   # TODO: can we do this in-place?
   for i = 1:k
      w1 = c[i] * vrot[i] - s[i] * vrot[i+1]
      # Here's a modest change that makes the code work in complex
      # arithmetic. Thanks to Howard Elman for this.
      # w2 = s(i)*vrot(i)+c(i)*vrot(i+1);
      w2 = s[i] * vrot[i] + c[i]' * vrot[i+1]
      vrot[i:i+1] = [w1, w2]
   end
   return vrot
end



"""
`function dgmres(f0, f, xc, errtol, kmax, [reorth, xinit]) -> x, error, total_iters`

GMRES linear equation solver for use in Newton-GMRES solver

### Input Parameters

* f0 : function at current point, f(xc)
* f : nonlinear function the format for f is  function fx = f(x)
              Note that for Newton-GMRES we incorporate any
              preconditioning into the function routine
* xc : current point
* errtol : relative residual reduction factor
* kmax : max number of iterations
* reorth : reorthogonalization method
      - 1 -- Brown/Hindmarsh condition (default)
      - 2 -- Never reorthogonalize (not recommended)
      - 3 -- Always reorthogonalize (not cheap!)
* xinit : initial iterate. xinit = zeros(f0) is the default. This
              is a reasonable choice unless restarted GMRES
              will be used as the linear solver.

### Output Parameters

* x : solution
* error : vector of residual norms for the history of the iteration
* total_iters : number of iterations
"""
function dgmres(f0, f, xc, errtol, kmax, reorth = 1, x = zeros(f0); h = 1e-7)
   # The right side of the linear equation for the step is -f0.
   b = - f0
   n = length(b)
   r = - dirder(xc, x, f, f0, h) - f0

   h = zeros(kmax, kmax)
   v = zeros(n, kmax)
   c = zeros(kmax+1)
   s = zeros(kmax+1)
   rho = norm(r)
   g = rho * eye(kmax+1, 1)
   errtol = errtol * norm(b)
   error = Float64[]

   # Test for termination on entry.
   push!(error, rho)
   total_iters = 0
   if (rho < errtol)
      # early termination
      return x, error, total_iters
   end

   v[:,1] = r / rho
   beta = rho
   k = 0

   # GMRES iteration
   while rho > errtol && k < kmax
      k = k+1;

      # Call directional derivative function.
      v[:, k+1] = dirder(xc, v[:,k], f, f0, h)
      normav = norm(v[:,k+1])

      # Modified Gram-Schmidt
      for j = 1:k
         h[j, k] = dot(v[:,j], v[:,k+1])
         v[:, k+1] = v[:,k+1] - h[j,k] * v[:,j]
      end
      h[k+1,k] = norm(v[:,k+1])
      normav2 = h[k+1,k]

      # Reorthogonalize?
      if  ((reorth == 1) && (normav + .001 * normav2 == normav)) || (reorth ==  3)
         for j = 1:k
            hr = dot(v[:,j], v[:,k+1])
            h[j,k] = h[j,k] + hr
            v[:,k+1] = v[:,k+1] - hr * v[:,j]
         end
         h[k+1, k] = norm(v[:,k+1])
      end

      # Watch out for happy breakdown.
      if h[k+1,k] != 0
         v[:,k+1] = v[:,k+1] / h[k+1,k]
      end

      # Form and store the information for the new Givens rotation.
      if k > 1
         h[1:k,k] = givapp(c[1:k-1], s[1:k-1], h[1:k,k], k-1)
      end

      # Don't divide by zero if solution has  been found.
      nu = norm(h[k:k+1,k])
      if nu != 0
         # c[k] = h[k,k]/nu
         c[k] = conj(h[k,k] / nu)
         s[k] = -h[k+1,k] / nu;
         h[k,k] = c[k] * h[k,k] - s[k] * h[k+1,k]
         h[k+1,k] = 0
         g[k:k+1] = givapp(c[k], s[k], g[k:k+1], 1)
      end

      # Update the residual norm.
      rho = abs(g[k+1])
      push!(error, rho)
   end # end of the main while loop

   # At this point either k > kmax or rho < errtol.
   # It's time to compute x and leave.
   y = h[1:k,1:k] \ g[1:k]
   total_iters = k
   x = x + v[1:n, 1:k] * y

   return x, error, total_iters, total_iters
end



export NK, blocklanczos


"""
`pushcol(V::Matrix, v::Vector) -> V`:
append the vector `v` as a column to `V`
"""
function pushcol(V::Matrix, v::Vector)
   rows, cols = size(V)
   V = V[:]
   append!(V, v)
   return reshape(V, rows, cols+1)
end


"""
`sorted_eig(A::Union{SymTridiagonal, Symmetric})`:
diagonalise A, and return sorted eigenvalues and analogously sorted eigenvectors
"""
function sorted_eig(A)
   D, Q = eig(A)
   I = sortperm(real(D))
   return D[I], Q[:, I]
end

nkdualnorm(P, f) = norm(f)  # TODO: revisit

"""
`orthogonalise(w, V::Matrix, P=I) -> w'`:
Assuming V is P-orthogonal, P-orthogonalise w to V and return; if
w is already in the span of V, then a random vector is orthogonalised
and returned instead. The returned vector is normalised.
"""
function orthogonalise(w, V::Matrix, P=I; maxiter = 10, ORTHTOL=1e-10)
   @assert size(V,2) < size(V,1)
   for n = 1:maxiter
      w /= norm(P, w)
      # make w orthogonal to all columns of V
      for i = 1:size(V,2)
         w -= dot(w, P, V[:,i]) * V[:,i]
      end
      # test that the resulting w is truly orthogonal:
      nrmw = norm(P, w)
      if nrmw < 1e-15 && n == maxiter; break; end
      if nrmw < 1e-15
         w = P \ (rand(length(w)) - 0.5)
         w /= norm(w, P)
      else
         w /= norm(P, w)
         if norm( w' * V, Inf ) < ORTHTOL
            return w
         end
      end
   end
   warn("""`orthogonalise` did not terminate succesfully, most likely
      the returned vector is not orthogonal to V""")
   return w
end


"""
`appendkrylov(V, AxV, Y, v, Hmul, P)`:
appends [V v], [AxV A*v], [Y P \ A*v]
"""
function appendkrylov(V, AxV, Y, v, Hmul, P)
   V   = pushcol(V,   v)
   AxV = pushcol(AxV, Hmul(v))
   Y   = pushcol(Y,   P \ AxV[:, end])
   return V, AxV, Y
end


"""
`immutable LanczosMatrix`: (todo: write documentation)
"""
immutable LanczosMatrix
   V::Matrix{Float64}
   D::Vector{Float64}
   P
end

Base.length(K::LanczosMatrix) = length(K.D)
Base.getindex(K::LanczosMatrix, i) = K.D[i], K.V[:,i]
# (Base.*)()  (TODO)
# collect
# \

"spectrum transformation for `dlanczos` for minimisation"
minim(D) = abs.(D)
"spectrum transformation for `dlanczos` for index-1 saddles"
index1(D) = [-abs(D[1]); abs.(D[2:end])]
"spectrum transformation for `dlanczos` for index-2 saddles"
index2(D) = [-abs.(D[1:2]); abs.(D[3:end])]



# TODO:
#   * look into re-orthogonalising
#   * look into getting rid of AxV
#   * this code does not make proper use of the lanczos recursion
#     at all; either switch to an arnoldi method, or clean it up

"""
dlanczos(f0, f, xc, errtol, kmax, transform; kwargs...) ->

Given E : ℝᴺ → ℝ with hessian H = ∇²E(x) and b ∈ ℝᴺ  `dlanczos` computes an
approximate solution `u` to the system
   g(H) u = b
where g(H) ∈ ℝᴺˣᴺ is a transformation of J = Q D Q'  of the form
   g(H) = Q g(D) Q'
Here, D = diag(λ₁, λ₂, …) is the diagonal matrix of ordered eigenvalues
(λ₁ being the smallest). More on the transformation g see below.

`dlanczos` first uses a (possibly preconditioned) (block-)lanczos method to compute
a reasonable approximation to H in the form H = P V T V' P. Then it diagonalises
T = Q D Q', replaces D with g(D) as above. This gives an approximation G to g(H).
The method terminates if `pinv(G) b` yields a solution of sufficiently high accuracy
(residual).

In the lanzcos iterations the matrix vector product H * u
is replaced with the finite-difference operation (∇E(xc + h u) - ∇E(xc))/h;
see also `dirder` and `dirderinf`.

## Required Parameters

* `f0` : ∇E(xc)
* `f` : function to evaluate ∇E
* `xc` : current point
* `b` : right-hand side in the equation
* `errtol` : max-norm residual tolerance for termination
* `kmax` : maximum number of lanzcos iterations
* `transform` : the g transformation, default: id

## KW parameters

* `P` : preconditioner (default: I;  minimally needs to define `*` and `\` )
* `V0` : initial subspace (default: P \ b)
* `eigatol`, `eigrtol`: tolerance on the first eigenvalue
* `debug`: show debug information (true/false)
* `h` : finite-difference parameter
* `dirder` : function to compute the directional derivative; see
         `dirderinf` for format
* `ORTHTOL` : orthogonality tolerance parameter
* `Hmul` : the finite-difference operation can in principle be replaced by
         and arbitrary operator specifying the matrix-vector product

## Returns

* x : approximate solution
* G : of type `LanczosMatrix`, to extract information about the computed operator
* success : true if termination criteria are satisfied, false if kmax is reached

## Further Notes

* Acknowledgement: This code is inspired by C T Kelley's Matlab Newton-Krylov solver `nsoli.m`;
see also `nsoli` and `dgmres` which are directly ported with permission from `nsoli.m`.
* Termination criterion: If g ≠ id then G * u - b is not in fact the residual
of the equation we are trying to solve but only an approximation to that residual.
Therefore the termination criterion is only approximately satisfied.
"""
function dlanczos( f0, f, xc, b, errtol, kmax, transform = identity;
                  P = I, V0 = P \ b,
                  eigatol = 1e-1, eigrtol = 1e-1,
                  debug = false, h = 1e-7,
                  dirder = dirderinf,
                  Hmul = z -> dirder(f, f0, xc, z, h),
                  ORTHTOL = 1e-12, nevals=1 )

   # make sure V0 is given in the correct format (Matrix)
   if isa(V0, Vector)
      V0 = reshape(V0, length(V0), 1)
   end

   # initialise some variables
   p = size(V0, 2)     # block size
   d = length(f0)      # problem dimension
   @assert kmax <= d
   numf = 0            # count f evaluations
   isnewton = false    # remember whether the output is a newton direction

   # initialise Krylov subspace and more
   V = zeros(d,0)      # store the Krylov basis
   AxV = zeros(d, 0)   # store A vⱼ
   Y = zeros(d, 0)     # store P \ A vⱼ

   # initialise Krylov subspace; TODO: detect if a vector is in the span of previous ones
   for j = 1:size(V0, 2)
      vj = orthogonalise(V0[:,j], V, P)
      V, AxV, Y = appendkrylov(V, AxV, Y, vj, Hmul, P)
   end

   # prepare for the Block-Lanczos loop
   j = 1
   x = zeros(d)
   vmin = zeros(d)
   λ = sort(eigvals(Symmetric(V'*AxV)))[1:nevals]   # TODO: possibly track more than one e-val
   λ_old = Vector{Float64}[]

   if debug
      @printf("      numf     λ    err_λ     |Ax-b|/|b| \n")
   end


   # start the block-lanczos loop; when we have kmax v-vectors we stop
   while size(V, 2) <= kmax

      # At this point we have V, AxV and Y = P \ AxV  available, so we can
      # now solve the projected linear system and eigenvalue problem
      #
      # there is probably an elegant way to assemble the projected
      # linear system, but for now we just do it brute-force:
      #    A vj = P V T V' P vj
      #    yj = P \ A qj = V T  V' P vj
      #    V' P Y = V' P V T V' P V = T
      n = size(V, 2)
      T = Symmetric(V' * AxV)       # T = V' * P * Y = V' * AV
      #  TODO: could replace above with with [dot(V[:,a], P, Y[:,b] for a=1:n, b=1:n]
      #        to avoid storing AxV altogether

      # make the specified transformation
      D, Q = sorted_eig(T)
      E = transform(D)
      # residual estimate for the old x    >>> TODO: should we switch to P^{-1}-norm?
      res = norm( P * (V * (Q * (E .* (Q' * (V' * (P * x)))))) - b )
      # new x and λ (remember the old)
      g = Q * (E .\ (Q' * (V' * b)))
      x, x_old = V * g, x
      push!(λ_old, λ)
      λ = D[1:neval]
      # if E == D then g(H) = H hence we can estimate the *actual* and *current* residual
      if (isnewton = (norm(E - D, Inf) < 1e-7) )
         res = norm( AxV * g - b )    # TODO: should we switch to P^{-1}-norm?
      end
      # check for termination
      err_λ = maximum(norm(λ - λ_old[i], Inf)
                      for i = max(1,length(λ)-p+1):length(λ))
      if debug
         @printf("      %d   %.2f   %.2e    %.2e \n", numf, λ, err_λ, res/norm(b))
      end

      # CHECK FOR TERMINATION
      if (res < errtol) && ((λ == E[1:neval]) || (err_λ < eigatol + eigrtol * abs(λ)))
         return x, LanczosMatrix(V * Q, E, P), true
      end

      # add the next Krylov vector
      if j <= size(V,2)   # we have Krylov vectors left to cycle through (99.9% of the time)
         w = orthogonalise(Y[:, j], V, P)
         V, AxV, Y = appendkrylov(V, AxV, Y, w/nrmw, Hmul, P)
         numf += 1
         j += 1
      else    # no vectors left, try to add a random vector (how can this happen???)
         warn("no vectors left; had to add a random vector")
         w = orthogonalise(P \ rand(d), V, P)
         V, AxV, Y = appendkrylov(V, AxV, Y, w/nrmw, Hmul, P)
         numf += 1
      end
   end
   # if we are here it means that kmax is reached, i.e. we terminate with
   # warning or error >>> TODO: return to how to best handle this?!
   # warn("`dcg_index1` did not converge within kmax = $(kmax) iterations")
   return x, λ, vmin, numf, isnewton
end



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
