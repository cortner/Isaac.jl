
# CTKSolvers

A repository for porting a subset of [Tim Kelley's](http://www4.ncsu.edu/~ctk/) Matlab code to Julia, with his kind permission, in particular the finite-difference Krylov methods and Newton-Krylov solvers.

Eventually this code will be integrated into [NLsolve.jl](https://github.com/JuliaNLSolvers/NLsolve.jl).

The initial translation is now complete and seems to work on all test problems. Please note that the current version is simply a direct translation of the Matlab code and as such has many inefficiencies that will (hopefully) be remedied over time.

## Installation

Since the package is not registered, use
```
Pkg.clone("git@github.com:cortner/CTKSolvers.jl.git")
```

## Usage Example

```
using CTKSolvers
f(x) = [x[1] + x[2] + x[1]^2, x[2] + x[1]*x[2]]
x, it_hist, ierr, x_hist = nsoli(rand(2), f)
```

For the full documentation see `?nsoli` in the REPL or IJulia



<!-- [![Build Status](https://travis-ci.org/cortner/CTKSolvers.jl.svg?branch=master)](https://travis-ci.org/cortner/CTKSolvers.jl)

[![Coverage Status](https://coveralls.io/repos/cortner/CTKSolvers.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/cortner/CTKSolvers.jl?branch=master)

[![codecov.io](http://codecov.io/github/cortner/CTKSolvers.jl/coverage.svg?branch=master)](http://codecov.io/github/cortner/CTKSolvers.jl?branch=master) -->
