using Test
using LinearAlgebra
using ppt2

function ispos(n::Int, m::Int; attempts=100000, atol=1e-6)
    C_phi = pncp_mat(n, m)

    # An arbitrary alternative Gram representative of the *same* generated form.
    Mλ = C_phi + sum(randn() * N for N in gram_freedom(n, m))

    for _ in 1:attempts
        xy = kron(randn(n), randn(m))

        @test xy' * C_phi * xy > -atol             # M₀ is (real) block positive
        @test xy' * Mλ * xy ≈ xy' * C_phi * xy     # Mλ represents the same polynomial
        @test xy' * Mλ * xy > -atol                # ... hence real block positive too
    end
end

pncp_mat(3, 3)

@testset "3x3 positive map" ispos(3, 3)
@testset "3x4 positive map" ispos(3, 4)
@testset "4x3 positive map" ispos(4, 3)
@testset "4x4 positive map" ispos(4, 4)
