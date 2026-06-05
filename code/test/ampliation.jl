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

# (Φ_A ⊗ I_m)(B) built block-wise from the Choi J(Φ_A): M_n → M_r — the
# `system=1` map, acting on the first subsystem.
function ampl_sys1_ref(A, B, n, m, r)
    out = zeros(eltype(A), r*m, r*m)
    for a in 1:n, c in 1:n
        Aac = A[(a-1)*r+1:a*r, (c-1)*r+1:c*r]      # Φ_A(E_ac)
        Bac = B[(a-1)*m+1:a*m, (c-1)*m+1:c*m]      # (a,c) block of B
        out += kron(Aac, Bac)                       # [r, m] ordering
    end
    return out
end

@testset "ampliation: system=1 (map on subsystem A)" begin
    # square: (Φ_A ⊗ I)(B) = SW · (I ⊗ Φ_A)(SW·B·SW) · SW
    for d in (2, 3, 4)
        A = randn(d^2, d^2); B = randn(d^2, d^2); S = ppt2.swap(d)
        @test isapprox(ampliation(A, B, d, d; system=1),
                       S * ampliation(A, S*B*S, d, d; system=2) * S, rtol=1e-10)
    end
    # rectangular: against the block-wise reference
    for (n, m, r) in ((2, 3, 4), (3, 2, 2), (2, 4, 3), (4, 2, 3))
        B = randn(n*m, n*m)        # operator on [n, m]
        A = randn(n*r, n*r)        # J(Φ_A): M_n → M_r
        C = ampliation(A, B, n, m; system=1)
        @test size(C) == (r*m, r*m)
        @test isapprox(C, ampl_sys1_ref(A, B, n, m, r), rtol=1e-9)
    end
    @test_throws ArgumentError ampliation(randn(4,4), randn(4,4), 2, 2; system=3)
    @test_throws DimensionMismatch ampliation(randn(5,5), randn(6,6), 2, 3; system=1)  # size(A) ∤ n
end

@testset "detect_ampliation: known positive-map detection" begin
    # Reduction map R(X) = tr(X)·I − X on M₂ (positive) detects the NPT Bell state.
    E(i, j) = (e = zeros(2, 2); e[i, j] = 1.0; e)
    R(X) = tr(X) * Matrix(I, 2, 2) - X
    C_R = sum(kron(E(i, j), R(E(i, j))) for i in 1:2, j in 1:2)   # Choi of R
    ψ = [1.0, 0, 0, 1] / sqrt(2); ρ = ψ * ψ'
    res = detect_ampliation(ρ, [C_R], 2, 2)                       # system=1 by default
    @test res.detected
    @test res.idx == 1
    @test res.value < 0
    # a separable state is not flagged
    a = randn(2); b = randn(2); σ = kron(a*a', b*b')
    @test !detect_ampliation(σ, [C_R], 2, 2).detected
end

@testset "detect_trace works on complex (Hermitian) states" begin
    ψ = [1.0, 0, 0, 1] / sqrt(2); P = ψ * ψ'
    W = Matrix(0.5I, 4, 4) - P                    # tr(W·P) = 0.5 - 1 = -0.5 < 0
    r = detect_trace(ComplexF64.(P), [W])         # complex τ must not break findmin
    @test r.detected
    @test r.value ≈ -0.5
    @test r.idx == 1
end
