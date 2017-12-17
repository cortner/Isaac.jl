

@testset "dlanczos" begin
println("Testing `dlanczos`")

# create A = -Δ , B = diag(V)   with V positive
d = 30
x = linspace(-1,1,d)
A = d^2 * SymTridiagonal(2*ones(d), -ones(d-1))
B = Diagonal(exp.(- x.^2) + 0.1)

# construct an index-1 system
μ = 30.0
H = A - μ * B
σ = sort(eigvals(full(H)))
@assert σ[1] < 0 && σ[2] > 0

b = ones(d) + 0.1 * sin.(x)

# define a linear function
fl = x -> H * x

# then try the same with a quadratic nonlinearity mixed in
#   (hessian remains the same!)
q = x -> [ [x[i]*x[i+1] for i = 1:d-1]; x[end] * x[1] ]
fq = x -> H * x + q(x)

# testing unpreconditioned lanczos
for (f, V0, msg) in [ (fl, reshape(b, d, 1), "Linear - Basic Lanczos"),
                      (fl, [b rand(d)], "Linear - Block Lanczos"),
                      (fq, reshape(b, d, 1), "Nonlinear - Basic Lanczos"),
                      (fq, [b rand(d)], "Nonlinear - Block Lanczos")
                     ]
   println("Testing: ", msg)
   x, G, numf, success = dlanczos(zeros(d), f, zeros(d), b, 1e-6, d; V0 = V0, debug = false)
   λ, v = G[1]
   @test success
   @test norm(H * x - b) < 1e-6
   @test norm(H * v - λ * v, Inf) < 1e-6
   @test abs(λ - σ[1]) < 1e-7
end

# next we add some preconditioning
P = 0.9 * sparse(A) + μ/2 * speye(d)
λP = minimum(eigvals(full(H), full(P)))
srand(12345)
vrand = rand(d)

for (V0, msg, numfo) in [ (reshape(P\b, d, 1), "Preconditioned Lanczos", 9),
                           ([P\b P\vrand], "Preconditioned Block-Lanczos", 14) ]
   println("Testing: ", msg)
   x, G, numf, success = dlanczos(zeros(d), fq, zeros(d), b, 1e-6, d;
                                P = P, debug = false, V0=V0)
   @test numf <= numfo
   @test success
   λ, v = G[1]
   @test norm(H * x - b) < 1e-6
   @test norm(H * v - λ * P * v, Inf) < 1e-6
   @test abs(λ - λP) < 1e-6
   if numf != numfo
      warn("numf was $numfo in original tests; now it is $numf")
   end
end

end  # @testset "blocklanczos"
