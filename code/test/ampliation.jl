using Test
using LinearAlgebra
using Ket
using ppt2

# `ampliation(A, B, n, m)` should equal the Choi matrix of the map composition,
# computed directly from the Choi–Jamiołkowski identity
#
#     J(Φ_A ∘ Φ_B) = Tr_{2,3}[ (J(Φ_A) ⊗ J(Φ_B)) (I ⊗ |e⟩⟨e| ⊗ I) ],
#     |e⟩ = Σ_i |i⟩ ⊗ |i⟩.
#
# (Salvaged from the old `playground` notebook, where it was checked by hand.)
function choi_composition(A, B, d)
    e = sum(kron(ket(i, d), ket(i, d)) for i in 1:d)
    return partial_trace(kron(A, B) * kron(I(d), kron(e * e', I(d))), [2, 3], [d, d, d, d])
end

# In the PPT² test a map is composed with itself, ampliation(M, M, d, d), so we
# check the self-composition (where the argument order is unambiguous).
@testset "ampliation equals Choi composition" begin
    for d in (2, 3, 4)
        M = randn(d^2, d^2)
        @test isapprox(ampliation(M, M, d, d), choi_composition(M, M, d), rtol=1e-10)
    end
end

# General (rectangular) link product for J(Φ_A ∘ Φ_B):
#   Φ_B: M_n → M_m  (B on systems 1,2),  Φ_A: M_m → M_r  (A on systems 3,4),
# linked over the shared m-system by |e⟩⟨e|, |e⟩ = Σ_i |i⟩⊗|i⟩ ∈ C^m ⊗ C^m.
function choi_composition(A, B, n, m, r)
    e = sum(kron(ket(i, m), ket(i, m)) for i in 1:m)
    M = kron(B, A) * kron(I(n), kron(e * e', I(r)))
    return partial_trace(M, [2, 3], [n, m, m, r])
end

@testset "ampliation: rectangular composition" begin
    for (n, m, r) in ((2, 3, 4), (3, 2, 2), (2, 4, 3), (4, 2, 3), (3, 3, 2))
        B = randn(n*m, n*m)        # J(Φ_B): M_n → M_m
        A = randn(m*r, m*r)        # J(Φ_A): M_m → M_r
        C = ampliation(A, B, n, m)
        @test size(C) == (n*r, n*r)
        @test isapprox(C, choi_composition(A, B, n, m, r), rtol=1e-9)
    end
end

@testset "ampliation: dimension checks" begin
    @test_throws DimensionMismatch ampliation(randn(6, 6), randn(6, 6), 2, 2)  # B not (n·m)²
    @test_throws DimensionMismatch ampliation(randn(5, 5), randn(6, 6), 2, 3)  # size(A) ∤ m
end
