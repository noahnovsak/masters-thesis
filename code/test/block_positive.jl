using Test
using LinearAlgebra
using Random
using ppt2

# `is_block_positive` samples product vectors and reports a violation of
# ⟨x⊗y|W|x⊗y⟩ ≥ 0. PSD matrices are block-positive over both fields. The point of
# these tests is the asymmetry introduced by `gram_freedom`: its directions vanish
# on *real* product vectors (so real block-positivity is preserved for any
# coefficient) but not on complex ones (so complex block-positivity — the genuine
# witness condition — can be destroyed). This is exactly why isolating the valid
# representatives `Mλ` is a real problem and `field=:real` cannot do it.
@testset "is_block_positive" begin
    rng = MersenneTwister(11)

    for (n, m) in ((2, 2), (2, 3), (3, 3))
        d = n * m

        # a PSD matrix is block-positive over ℝ and ℂ
        P = rand_psd(n, m; rng = rng)
        @test is_block_positive(P, n, m; field = :real,    trials = 2000, rng = rng)
        @test is_block_positive(P, n, m; field = :complex, trials = 2000, rng = rng)

        # a negative-definite matrix is rejected immediately
        @test !is_block_positive(-Matrix(I, d, d), n, m; field = :complex, trials = 10, rng = rng)

        # Mλ = P + t·N stays real block-positive for any t (N vanishes on real
        # products) but loses complex block-positivity once t is large.
        N = gram_freedom(n, m)[1]                       # the (i,j,k,l) = (1,1,2,2) relation
        t = 1e3 * (opnorm(P) + 1)
        W = P + t * N

        @test is_block_positive(W, n, m; field = :real, trials = 2000, rng = rng)

        # explicit complex product vector violating ⟨z|W|z⟩ ≥ 0. For the (1,1,2,2)
        # relation, ⟨z|N|z⟩ = -4·Im(x̄₁x₂)·Im(ȳ₁y₂); choosing x₁=y₁=1, x₂=y₂=i gives -4.
        x = zeros(ComplexF64, n); x[1] = 1; x[2] = im
        y = zeros(ComplexF64, m); y[1] = 1; y[2] = im
        z = kron(x, y)
        @test real(z' * N * z) ≈ -4
        @test real(z' * W * z) < 0
        @test !is_block_positive(W, n, m; field = :complex, trials = 5000, rng = rng)
    end
end
