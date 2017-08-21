
using SaddleSearch
using SaddleSearch.TestSets
using Plots, Dierckx, LaTeXStrings


function contourPlot{T}(x::Vector{T})
    pyplot(leg=false, ticks=nothing)

    dX = linspace(-1.2, 1.0, 100) |> collect
    dY = linspace(-0.2, 1.7, 100) |> collect
    X = repmat(dX', length(dY), 1)
    Y = repmat(dY, 1, length(dX))
    E(dX,dY) = begin
        Energy([dX; dY])
    end
    Z = map(E, X, Y)
    p = Plots.contour(dX, dY, E, fill=true, color=:heat, grid=false)

    xrefined=[]
    ds = [norm(x[i+1]-x[i]) for i=1:length(x)-1]
    s = [0; [sum(ds[1:i]) for i in 1:length(ds)]]; s /= s[end]; s[end] = 1.
    S = [Spline1D(s, [x[j][i] for j=1:length(s)],
        w = ones(length(x)), k = 3, bc = "error") for i=1:length(x[1])]
    xrefined = cat(2,[[S[i](s) for i in 1:length(S)] for s in linspace(0., 1., 120) ]...)
    p = Plots.plot!(xrefined[1,:], xrefined[2,:])

    xmat = cat(2, x...)
    p = Plots.plot!(xmat[1,:], xmat[2,:], marker=(4.0,:c), linewidth=0.)

    return p
end

function convergencePlot(log, label, colour, style)
  evaluations = log.D[:numdE]
  res = log.D[:maxres]

  Plots.plot( evaluations, res, yscale=:log10, label=label,
            linecolor=colour, linestyle=style, linewidth=2,
            xlabel="force evaluations",
            ylabel=(L"$\mathrm{sup}_k||\nabla E^\perp||$") )
end

V = MullerPotential()
Energy, Gradient = objective(V)
d = 2

# ==================== STANDARD FIXED STEP SOLVE ==================
# SET-UP PARAMETERS
α = 0.0001 # fixed step size(very slow, but converging with this step size)
tol = 1e-4 # force residual tolerance limit (termination criterion)
maxnit = 500 # maximum number of iterations
verbose = 1
precon_cond = false

path = SaddleSearch.StringMethod(α, tol, maxnit, I, (P, x) -> P, verbose, precon_cond)

# INITIAL LINEAR PATH
N = 25 # number of images
x0 = [-0.5; 1.5]; x1 = [0.7; .0] # local minima
x = [(1-s)*x0 + s*x1 for s in linspace(.0, 1., N)] # the path
t = [((x1-x0)/norm(x1-x0)) for s in linspace(.0, 1., N)]; # tangent along path

PATHx, PATHlog = SaddleSearch.run!(path, Energy, Gradient, x, t)

contourPlot(PATHx)

# convergencePlot(PATHlog, ("α = $α"), :blue, :solid)



# ==================== NEWTON SOLVE ==================

function get_lam_x(lamx)
    N = (1 + length(lamx)) ÷ (1 + d)
    λ = lamx[1:N-1]
    x = reshape(lamx[N:end], d, N)
    return N, λ, x
end

function fstring1(lamx)
    N, λ, x = get_lam_x(lamx)
    F = zeros(d, N)
    for n = 1:N
        F[:, n] = Gradient(x[:,n])
    end
    for n = 2:N-1
        t = x[:,n+1]-x[:,n-1]
        F[:,n] -= λ[n] * t / norm(t)
    end
    C = [ norm(x[:,i]-x[:,i-1]) for i = 2:N ]
    return [F[:]; C - λ[1]]
end

function project(lamx)
    # here we should probably redistribute
    N, λ, x = get_lam_x(lamx)
    λ[1] = mean(  norm(x[:,i]-x[:,i-1]) for i = 2:N  )
    for n = 2:N-1
        t = x[:,n+1]-x[:,n-1]
        dE = Gradient(x[:,n])
        λ[n] = dot(t, dE)
    end
    return [λ; x[:]]
end

using Isaac

z = zeros(d, N)
for n = 1:N
    z[:,n] = x[n]
end
lamx0 = project( [zeros(N-1); z[:]] )

lamx0, nde = Isaac.nsolistab(fstring1, lamx0; maxstep = 0.1, verbose = 3)







##### BACKUP
# spline(s, y) = Spline1D(s, y, w = ones(length(y)), k=3, bc="error")
#
# function reparametrise(x)
#     N = length(x); d = length(x[1])
#     ds = [norm(x[i+1]-x[i]) for i=1:N-1]
#     s = [0.0; cumsum(ds)]
#     s /= s[end]; s[end] = 1.0
#     ŝ = linspace(0., 1., N)
#     S = [ spline(s, [x[j][i] for j=1:N]) for i=1:d ]
#     x = [ [S[i](r) for i = 1:N] for r in ŝ ]
#     t = [ [derivative(S[i], r) for i in 1:N] for r in ŝ ]
#     return x, t
# end
#
# function fstring(z)
#     # convert to list of nodes
#     z = reshape(z, d, length(x) ÷ d)
#     x = [ z[:,j] for j = 1:size(z,2) ]
#     # reparametrise and compute tangents
#     x, t = reparametrise(x)
#     # compute the force on the string
#
#     dE0 = [Gradient(x[i]) for i=1:length(x)]
#     dE0⟂ = [dE0[i] - dot(dE0[i],t[i])*t[i] for i = 1:length(x)]
# end
