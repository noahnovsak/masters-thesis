using Test
using LinearAlgebra
using ppt2

# `swap` and `symmetric_projector` are internal (not exported, since `Ket`
# already exports `symmetric_projector`), so reach them as `ppt2.<name>`.

@testset "swap operator" begin
    for d in (2, 3, 4)
        S = ppt2.swap(d)
        @test size(S) == (d^2, d^2)
        @test S == S'                       # symmetric
        @test S * S ≈ I                     # involution
        for _ in 1:5                        # swaps the factors of a product vector
            a = randn(d); b = randn(d)
            @test S * kron(a, b) ≈ kron(b, a)
        end
    end
end

@testset "(anti)symmetric projectors" begin
    for d in (2, 3, 4)
        P = ppt2.symmetric_projector(d)
        Q = antisymmetric_projector(d)
        @test P ≈ P'                        # Hermitian
        @test Q ≈ Q'
        @test P * P ≈ P                     # idempotent
        @test Q * Q ≈ Q
        @test P + Q ≈ I                     # complementary
        @test P * Q ≈ zeros(d^2, d^2)       # orthogonal ranges
        @test rank(P) == d * (d + 1) ÷ 2    # dim of symmetric subspace
        @test rank(Q) == d * (d - 1) ÷ 2    # dim of antisymmetric subspace
        S = ppt2.swap(d)                    # images are the ±1 eigenspaces of swap
        @test S * P ≈ P
        @test S * Q ≈ -Q
    end
end
