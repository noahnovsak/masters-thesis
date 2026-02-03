using LinearAlgebra
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using ppt2

include("./candidate.jl")
using .candidate

attempts = 20
entangled = true
d = 4
a = 0.5
t = 0.1
C = nothing

if entangled
    if d == 3
        C = C3()
    elseif d == 4
        C = C4(a, t)
    end
else
    C = ones(d^2, d^2)
end

println("Testing candidate: d = $d, $(entangled ? "entangled" : "separable")")

for i in 1:attempts
    del, v, V = gen_pncp(d, d)
    poly = del * v + 10 * vec(V * V')

    ent = reshape(poly, d^2, d^2) * C

    v_min = minimum(real(eigvals(ent)))

    println("Iteration $i: min eigenvalue = $v_min")
    if v_min < -1e-6
        break
    end
end
