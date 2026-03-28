using Test
using LinearAlgebra
using ppt2

function rand_sep(n::Int, m::Int; n_terms::Int=2)
    d = n * m
    rho = zeros(d, d)

    for _ in 1:n_terms
        psiA = randn(n)
        psiB = randn(m)

        rhoA = psiA * psiA'
        rhoB = psiB * psiB'

        rho += kron(rhoA, rhoB)
    end

    return rho
end

function ispos(n::Int, m::Int; attempts=100000, atol=1e-6)
    C_phi = pncp_mat(n, m)

    for _ in 1:attempts
        x = randn(n)
        y = randn(m)
        xy = kron(x, y)

        @test xy' * C_phi * xy > -atol
    end
end

function istrpos(n::Int, m::Int; attempts=100000, atol=1e-6)
    C_phi = pncp_mat(n, m)

    for _ in 1:attempts
        xy = rand_sep(n, m)

        @test tr(C_phi * xy) > -atol
    end
end

pncp_mat(3, 3)

@testset "3x3 positive" ispos(3, 3)
@testset "3x3 witness" istrpos(3, 3)
@testset "3x4 positive" ispos(3, 4)
@testset "3x4 witness" istrpos(3, 4)
@testset "4x3 positive" ispos(4, 3)
@testset "4x3 witness" istrpos(4, 3)
@testset "4x4 positive" ispos(4, 4)
@testset "4x4 witness" istrpos(4, 4)
