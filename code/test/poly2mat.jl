using Test
using LinearAlgebra
using Random
using ppt2

# `poly2mat(c, n, m)` builds the Gram matrix M (in the z = X⊗Y monomial basis)
# of the quartic form  c ⋅ ((X⊗Y) ⊗ (X⊗Y)).  Its defining property is therefore
#
#     zᵀ M z  ==  c ⋅ (z ⊗ z),     z = x ⊗ y,
#
# for every numeric x ∈ ℝⁿ, y ∈ ℝᵐ.  The matrix must also be symmetric.
@testset "poly2mat conversion" begin
    rng = MersenneTwister(42)
    for (n, m) in ((2, 2), (2, 3), (3, 3), (3, 4))
        d = n * m
        c = randn(rng, d^2)
        M = ppt2.poly2mat(c, n, m)
        @test size(M) == (d, d)
        @test M ≈ M'
        for _ in 1:10
            x = randn(rng, n); y = randn(rng, m)
            z = kron(x, y)
            @test z' * M * z ≈ c ⋅ kron(z, z)
        end
    end
end
