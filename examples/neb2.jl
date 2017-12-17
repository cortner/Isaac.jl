

path = @__FILE__()[1:end-7]

using SaddleSearch, SaddleSearch.TestSets, Isaac, JLD

V = MullerPotential()
Energy, Gradient = objective(V)
d = 2


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

function fneb(x, k, gradE)
    X = dof2vecs(x)
    N = length(X)
    dE = [gradE(X[i]) for i=1:N]
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

Dfneb(x, k, gradE) = ForwardDiff.jacobian(z -> fneb(z, k, gradE), x)

# setup NEB with spring constant 1 (scaled against number of nodes)
k = 1.0
Fneb(x) = fneb(x, k, Gradient)
Fneb0(x) = fneb(x, 0.0, Gradient)
DFneb(x) = Dfneb(x, k, Gradient)
Pneb(x, μ) = Dfneb(x, k, x->μ*x)

# load a very good starting guess so we can be sure (hope) to be in the
# region of attraction of Newton's method
z = JLD.load(path * "near.jld", "z")
N = length(z) ÷ d

# z = reshape(z, d, N)
# z = z[:, 1:2:end]
# z = z[:]
# N = length(z) ÷ d

println("Trying a few Newton steps")
x = copy(z)
@show norm(Fneb(x), Inf)
for n = 1:7
    x -= DFneb(x) \ Fneb(x)
    @show norm(Fneb(x), Inf)
end

println("Trying NSOLI")
x = copy(z)
z_nsoli, nde, ierr, xhist = Isaac.nsoli(x, Fneb, atol=1e-5, rtol=1e-5)
numF = round(Int, nde[end,2])
@show norm(Fneb(z_nsoli), Inf)
@show numF


println("Trying NSOLI-STAB")
x = copy(z)
z_nsoli, numF = Isaac.nsolistab(Fneb, x, tol=1e-5)
@show norm(Fneb(z_nsoli), Inf)
@show numF


function Pdiag(x)
    N = length(x) ÷ d
    P = zeros(length(x), length(x))
    for n = 1:N
        In = ((n-1) * 2) .+ (1:2)
        h = ForwardDiff.jacobian(Gradient, x[In])
        λ, V = eig(h)
        P[In, In] = V * diagm(abs.(λ)) * V'
    end
    return P
end


println("Trying P-NSOLI-STAB")
x = copy(z)
P = y -> Pdiag(y) + Pneb(y, 0.0)
z_nsoli, numF = Isaac.nsolistab(Fneb, x, tol=1e-5,
            P = P(x), precon_prep = (P_, x) -> P(x))
@show norm(Fneb(z_nsoli), Inf)
@show numF

println("Trying Preconditioned NSOLI")
x = copy(z)
pFneb = y ->  P(y) \ Fneb(y)
z_nsoli, nde, ierr, xhist = Isaac.nsoli(x, pFneb, atol=1e-7, rtol=1e-7)
numF = round(Int, nde[end,2])
@show norm(Fneb(z_nsoli), Inf)
@show numF



# using Plots
#
# J = DFneb(z)
# p = P(z)
# λ, V = eig(p \ J)
# plot(real(λ), imag(λ), lw = 0, m = :o)
#
# minimum(real.(λ))



# bad starting guess

d = 2
z = JLD.load(path * "far.jld", "z")

x = copy(z)
z_nsoli, numF = Isaac.nsolistab(Fneb, x, tol=1e-5, maxstep = 0.1, verbose = 3)
@show norm(Fneb(z_nsoli), Inf)
@show numF
