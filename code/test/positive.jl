using Test
using ppt2

function ispos(n::Int, m::Int, attempts::Int, atol::Float64)
    C_phi = pncp_form(n, m)

    for _ in 1:attempts
        x = randn(n)
        y = randn(m)
        xy = kron(x, y)

        @test xy' * C_phi * xy > -atol
    end
end

gen_pncp(3, 3)

@testset "generate 3x3 positive form" ispos(3, 3, 1000000, 0.0)
@testset "generate 3x4 positive form" ispos(3, 4, 1000000, 0.0)
@testset "generate 4x3 positive form" ispos(4, 3, 1000000, 0.0)
@testset "generate 4x4 positive form" ispos(4, 4, 100000, 0.0)
