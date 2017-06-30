
module CTKSolvers

export nsoli, dgmres

# linsearch  method
include("linesearch.jl")

# specialised Krylov subspace methods using finite-differences
include("fdkryloc.jl")

# the actual Newton-Krylov solver
include("nsoli.jl")


end
