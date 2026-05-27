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
