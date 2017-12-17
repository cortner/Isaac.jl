

path = @__FILE__()[1:end-6]

using SaddleSearch
using SaddleSearch.TestSets
using Plots, Dierckx, LaTeXStrings
using Isaac
using JLD

function contourPlot{T}(x::Vector{T}; v=nothing, text = nothing)
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
    xmat = cat(2, x...)
    p = Plots.plot!(xmat[1,:], xmat[2,:], marker=(4.0,:c), linewidth=0.)

    if v != nothing
        Plots.quiver!(xmat[1,:], xmat[2,:], quiver = (v[1,:], v[2,:]))
    end

    if text != nothing
        p = annotate!(-0.2, -0.7, text)
    end

    return p
end

V = MullerPotential()
Energy, Gradient = objective(V)
d = 2


z = JLD.load(path * "near.jld", "z")
contourPlot(dof2vecs(z))

# ================ NEB =================

vecsign(x) = x / norm(x)

function dof2vecs(x)
    N = length(x) ÷ d
    x = reshape(x, d, N)
    return [ x[:, n] for n = 1:N ]
end

function vecs2dof(X)
    x = X[1]
    for n = 2:length(X)
        append!(x, X[n])
    end
    return x
end

function fneb(x, k)
    X = dof2vecs(x)
    N = length(X)
    dE = [Gradient(X[i]) for i=1:N]
    # central finite differences
    dxds = [ vecsign(X[i+1]-X[i-1]) for i=2:N-1 ]
    dxds = [ [zeros(d)]; dxds; [zeros(d)] ]
    # elastic
    Fk = k*N^2*[dot(X[i+1] - 2*X[i] + X[i-1], dxds[i]) * dxds[i] for i=2:N-1]
    Fk = [[zeros(d)]; Fk; [zeros(d)]]
    # nudge
    dE0⟂ = [dE[i] - dot(dE[i],dxds[i])*dxds[i] for i = 1:N]
    # return total force
    return vecs2dof(dE0⟂ - Fk)
end

Dfneb(x, k) = ForwardDiff.jacobian(z -> fneb(z,k), x)


k = 1.0
Fneb(x) = fneb(x,k)
DFneb(x) = Dfneb(x, k)

z = JLD.load(path * "near.jld", "z")
N = length(z) ÷ d

# @show norm(Fneb(z), Inf)
# contourPlot(dof2vecs(z))


# J = DFneb(z)
# Λ, V = eig(J)
# anim = @animate for a = 1:length(Λ)
#     vec = reshape(V[:,a] |> real, 2, N)
#     p = contourPlot(dof2vecs(z), v = vec, text = "λ$a=$(Λ[a])")
#     p = Plots.plot!(title = "eigen-spectrum of Jacobian")
# end
# Plots.gif(anim, "neb_Jevecs.gif", fps=1)

# plot(real(Λ), imag(Λ), lw=0, marker=:o)
# plot!([0.0], [0.0], lw=0, marker = (:hex, :red, 10))
# extrema(abs(Λ))

# try a few newton steps with pseudo-inverse
x = copy(z)
@show norm(Fneb(x), Inf)
for n = 1:7
    x -= DFneb(x) \ Fneb(x)
    @show norm(Fneb(x), Inf)
end

contourPlot(dof2vecs(x))

x = copy(z)
z_nsoli, nde, ierr, xhist = Isaac.nsoli(x, Fneb)
numF = size(nde, 1)
@show norm(Fneb(z_nsoli), Inf)
@show numF, ierr
