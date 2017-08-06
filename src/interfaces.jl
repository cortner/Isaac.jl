using Parameters

# """
# Newton-Krylov based saddle search method
# """
# @with_kw type NK
#    tol::Float64 = 1e-5
#    maxnumdE::Int = 1000
#    len::Float64 = 1e-7
#    precon = I
#    precon_prep! = (P, x) -> P
#    verbose::Int = 1
#    krylovinit::Symbol = :resrot  # allow res, rand, rot, resrot
#    maxstep::Float64 = Inf
#    eigatol::Float64 = 1e-1
#    eigrtol::Float64 = 1e-1
# end
