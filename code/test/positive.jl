using Test
using ppt2

function ispos(n::Int, m::Int, attempts::Int, atol::Float64)
    del, v, V = gen_pncp(n, m)

    phi = del * v + 10 * vec(V * V')

    C_phi = vec(poly2mat(phi, n, m))

    for _ in 1:attempts
        x = randn(n)
        y = randn(m)
        xy = kron(x, y)

        if C_phi' * kron(xy, xy) < -atol
            return false
        end
    end

    return true
end

@testset "Sanity check: test positivity at random points" begin
    @test ispos(3, 3, 10000, 1e-6)
    @test ispos(3, 4, 10000, 1e-6)
    @test ispos(4, 3, 10000, 1e-6)
    @test ispos(4, 4, 10000, 1e-6)
    @test ispos(5, 5, 1000, 1e-6)
end
