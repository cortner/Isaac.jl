
module Isaac

include("testsets.jl")

# linsearch  method
include("linesearch.jl")

# specialised Krylov subspace methods using finite-differences
include("fdkrylov.jl")

# the Newton-Krylov solver based on C T Kelley's nsoli.m
include("nsoli.jl")

# the main transformed Newton solver
include("solvermain.jl")

# # interfaces to Optim.jl and NLsolve.jl
# include("interfaces.jl")


end
