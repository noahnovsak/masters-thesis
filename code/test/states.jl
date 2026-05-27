using Test
using LinearAlgebra
using Random
using ppt2

issym(A; atol=1e-9) = isapprox(A, A'; atol=atol)
ispsd(A; tol=1e-8) = eigmin(Hermitian(Matrix(A))) ≥ -tol

@testset "rand_psd" begin
    rng = MersenneTwister(1)
    for (n, m) in ((2, 2), (2, 3), (3, 3))
        d = n * m
        M = rand_psd(n, m; rng=rng)
        @test size(M) == (d, d)
        @test issym(M)
        @test ispsd(M)
        @test rank(M) == d                       # full rank by default
        for r in 1:d
            Mr = rand_psd(n, m; r=r, rng=rng)
            @test issym(Mr)
            @test ispsd(Mr)
            @test rank(Mr) == r                  # rank-r construction
        end
    end
end

@testset "rand_sep is PSD and PPT" begin
    rng = MersenneTwister(2)
    for (n, m) in ((2, 2), (2, 3), (3, 2), (3, 3))
        ρ = rand_sep(n, m; n_terms=3, rng=rng)
        @test size(ρ) == (n * m, n * m)
        @test issym(ρ)
        @test ispsd(ρ)
        @test is_ppt(ρ, n, m)                    # separable ⟹ PPT (kron(A,B) layout)
    end
end

@testset "rand_ppt is PSD and PPT" begin
    rng = MersenneTwister(3)
    for (n, m) in ((2, 2), (2, 3), (3, 3), (3, 4))
        ρ = rand_ppt(n, m; rng=rng)
        @test size(ρ) == (n * m, n * m)
        @test issym(ρ)
        @test ispsd(ρ)
        # rand_ppt is PPT under the same [n, m] = kron(A_n, B_m) bipartition as
        # rand_sep / is_ppt / Ket. The rectangular case is the discriminating one:
        # it fails if the block structure is built for the [m, n] ordering instead.
        @test is_ppt(ρ, n, m)
    end
end

@testset "rand_ppt is reproducible" begin
    for (n, m) in ((2, 2), (3, 4))
        @test rand_ppt(n, m; rng=MersenneTwister(7)) == rand_ppt(n, m; rng=MersenneTwister(7))
    end
end

@testset "is_ppt detects NPT entanglement" begin
    # |Φ⁺⟩ = (|00⟩ + |11⟩)/√2 on C²⊗C²: its partial transpose has a negative
    # eigenvalue, so the maximally entangled state is NOT PPT.
    ψ = [1.0, 0, 0, 1] / sqrt(2)
    @test !is_ppt(ψ * ψ', 2, 2)
    # a pure product state is PPT
    a = randn(2); b = randn(2)
    @test is_ppt(kron(a * a', b * b'), 2, 2)
end
