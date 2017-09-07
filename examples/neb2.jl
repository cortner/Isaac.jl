

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

# setup NEB with spring constant 1 (scaled against number of nodes)
k = 1.0
Fneb(x) = fneb(x,k)
DFneb(x) = Dfneb(x, k)

# load a very good starting guess so we can be sure (hope) to be in the
# region of attraction of Newton's method
z = JLD.load(path * "near.jld", "z")
N = length(z) ÷ d

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
