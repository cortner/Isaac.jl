
# ModulatedNewtonMethods.jl

This package provides a small set of Jacobian-Free Newton-Krylov solvers,
specifically:
* `nsoli` : a port of a subset of [Tim Kelley's](http://www4.ncsu.edu/~ctk/) `nsoli.m` Matlab code to Julia (with his kind permission). Any bugs or errors are of course my own.
* `nsolimod` : A *Modulated Jacobian-Free Newton-Krylov Solver*; instead of comuting general critical points of an energy it computes only critical points of a specific spectrum signature. For example only minima, or only index-1 saddles.

While `nsoli` (or at least `nsoli.m`) is a robust and well-tested code, `nsolimod` is still experimental, and in particular comes with no theoretical convergence guarantee of any kind. However, initial tests indicate that it is both robust and performs well in particular when the energy landscape is not highly nonlinear (but could still be highly ill-conditioned).

## Getting Started

Since the package is not registered, use
```
Pkg.clone("git@github.com:cortner/CTKSolvers.jl.git")
```

There is no documentation, but the inline documentation is fairly extensive;
start with
```
?nsoli
?nsolimod
```
in the REPL or IJulia

<!-- Eventually this code will be integrated into [NLsolve.jl](https://github.com/JuliaNLSolvers/NLsolve.jl). -->


## Examples

```
using CTKSolvers
f(x) = [x[1] + x[2] + x[1]^2, x[2] + x[1]*x[2]]
x, it_hist, ierr, x_hist = nsoli(rand(2), f)
```

See tests for more examples.

## TODO

* Performance tuning
* current implementation of `dlanzcos` is very naive and doesn't exploit the usual structures in Lanczos iterations; this should eventually be fixed
* integrate with `Optim.jl`, `NLsolve.jl`, `NLSolversBase.jl`
* generalise codes to admit actual hessian-vector products and full hessian inversion
* add more examples


<!-- [![Build Status](https://travis-ci.org/cortner/CTKSolvers.jl.svg?branch=master)](https://travis-ci.org/cortner/CTKSolvers.jl)

[![Coverage Status](https://coveralls.io/repos/cortner/CTKSolvers.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/cortner/CTKSolvers.jl?branch=master)

[![codecov.io](http://codecov.io/github/cortner/CTKSolvers.jl/coverage.svg?branch=master)](http://codecov.io/github/cortner/CTKSolvers.jl?branch=master) -->
