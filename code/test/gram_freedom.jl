using Test
using LinearAlgebra
using Random
using ppt2

# The Gram matrix of a biquadratic form is not unique: `poly2mat` returns one
# representative `M₀`, and adding any element of the Segre-vanishing space `L`
# (spanned by `gram_freedom`) yields another matrix `Mλ` that represents the
# *same* polynomial. These tests check that `M₀` and arbitrary `Mλ` agree as
# forms on real product vectors (their defining property), and that the basis
# of `L` has the expected size, is independent, and indeed vanishes on products.
@testset "gram_freedom: representation freedom" begin
    rng = MersenneTwister(7)
    for (n, m) in ((2, 2), (2, 3), (3, 3), (3, 4))
        d = n * m
        B = gram_freedom(n, m)

        # dim L = C(n,2)·C(m,2), and the basis is linearly independent
        @test length(B) == binomial(n, 2) * binomial(m, 2)
        @test rank(reduce(hcat, vec.(B))) == length(B)

        for N in B
            @test N ≈ N'                                  # symmetric
            for _ in 1:10                                 # vanishes on real products
                z = kron(randn(rng, n), randn(rng, m))
                @test abs(z' * N * z) < 1e-10
            end
        end

        # M₀ from a single (random) polynomial, plus a random representative Mλ
        c  = randn(rng, d^2)
        M0 = ppt2.poly2mat(c, n, m)
        Mλ = M0 + sum(randn(rng) * N for N in B)
        @test Mλ ≈ Mλ'

        # Mλ reproduces the same polynomial as M₀ (equivalently as the coeffs c)
        for _ in 1:50
            z = kron(randn(rng, n), randn(rng, m))
            @test z' * Mλ * z ≈ z' * M0 * z
            @test z' * Mλ * z ≈ c ⋅ kron(z, z)
        end
    end
end
