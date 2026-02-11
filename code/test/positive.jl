using Test
using ppt2

function ispos(d::Int, attempts::Int, atol::Float64)
    del, v, V = gen_pncp(d, d)

    form = del * v + 10 * vec(V * V')

    M = vec(poly2mat(form, d, d))

    for _ in 1:attempts
        x = randn(d)
        y = randn(d)
        xy = kron(x, y)

        if M' * kron(xy, xy) < -atol
            return false
        end
    end

    return true
end

# precompile
ispos(3, 1, 1e-6)

@testset "Sanity check: test positivity at random points" begin
    @test ispos(3, 10000, 1e-6)
    @test ispos(4, 10000, 1e-6)
end
