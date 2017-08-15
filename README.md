
# ModulatedNewtonMethods.jl

<!-- [![Build Status](https://travis-ci.org/cortner/CTKSolvers.jl.svg?branch=master)](https://travis-ci.org/cortner/CTKSolvers.jl)

[![Coverage Status](https://coveralls.io/repos/cortner/CTKSolvers.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/cortner/CTKSolvers.jl?branch=master)

[![codecov.io](http://codecov.io/github/cortner/CTKSolvers.jl/coverage.svg?branch=master)](http://codecov.io/github/cortner/CTKSolvers.jl?branch=master) -->


This package provides a small set of Jacobian-Free Newton-Krylov solvers,
specifically:
* `nsoli` : a port of a subset of [Tim Kelley's](http://www4.ncsu.edu/~ctk/) `nsoli.m` [Matlab code](http://www4.ncsu.edu/~ctk/newton/SOLVERS/nsoli.m) to Julia (with his permission), which is a well-tested Jacobian-Free Newton Krylov Solver. Any bugs or errors are of course my own.
* `nsolimod` : A *Modulated Jacobian-Free Newton-Krylov Solver* : instead of computing general critical points of an energy or roots of a nonlinear system this solver computes only critical points of a specific spectrum signature. For example only minima, or only index-1 saddles.

## Getting Started

Since the package is not registered, use
```
Pkg.clone("git@github.com:cortner/CTKSolvers.jl.git")
```

There is no documentation, but the inline documentation is fairly extensive;
start with
```
using CTKSolvers
?nsoli
?nsolimod
```
in the REPL or IJulia


## Examples

See the `tests` directory for more examples.

### Nonlinear Solver `nsoli`

```
using CTKSolvers
f(x) = [x[1] + x[2] + x[1]^2, x[2] + x[1]*x[2]]
x, it_hist, ierr, x_hist = nsoli(rand(2), f)
```

### Minimisation and saddle-search with `nsolimod`

```
using CTKSolvers, ForwardDiff
E(x) = (x[1]^2 - 1)^2 + (x[2]-x[1]^2)^2
dE(x) = ForwardDiff.gradient(E, x)
# minimisation (index-0 saddle search)
xmin, nmin = nsolimod(dE, [0.6,0.1], 0)
# index-1 saddle-search
xsad, nsad = nsolimod(dE, [0.4,-0.1], 1)
```


## When should I use `ModulatedNewtonMethods`?

For most users looking for a robust and well-tested optimiser or nonlinear solver [`Optim.jl`](https://github.com/JuliaNLSolvers/Optim.jl) and [`NLsolve.jl`](https://github.com/JuliaNLSolvers/NLsolve.jl) are probably better choices. That said, the code `nsoli` (or at least its parent `nsoli.m`) is a robust and well-tested code, and always worth comparing against `NLsolve.jl` to decide which will perform better on a specific problem.

The second code, `nsolimod` is still very much experimental, and in particular comes with no theoretical convergence guarantee of any kind. However, initial tests indicate that it is both robust and performant, in particular when the energy landscape is not highly nonlinear but possibly ill-conditioned.

`nsolimod` is written with expensive objectives in mind where each gradient (or function) evaluation far outweighs the cost of the optimisation boilerplate and of the linear algebra. Over time, I hope to optimise the code so that it becomes competitive for cheap objectives.

The main use case for `nsolimod` at the moment is saddle-search. On my tests systems `nsolimod` outperforms all algorithms I have tested against, both in terms of performance and robustness.  (though I have not yet exhausted the space of all saddle search methods)

## TODO

* Rename and register the repository
* generalise `nsolimod` to nonlinear systems
* improve linesearch, especially for minimisation
* Performance tuning
* current implementation of `dlanzcos` is very naive as it does not exploit the usual structures in Lanczos iterations; this should eventually be fixed or maybe just move to Arnoldi
* integrate with `Optim.jl`, `NLsolve.jl`, `NLSolversBase.jl`
* generalise codes to admit actual hessian-vector products and full hessian inversion, depending on whether the objective is `OnceDifferentiable` or `TwiceDifferentiable`
* add more examples
