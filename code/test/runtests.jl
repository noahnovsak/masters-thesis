using Test
using JuMP
using MosekTools

# SDP-backed tests need a working Mosek license; skip them gracefully otherwise.
function mosek_licensed()
    try
        model = Model(Mosek.Optimizer)
        set_silent(model)
        @variable(model, x >= 1)
        @objective(model, Min, x)
        optimize!(model)
        return termination_status(model) == OPTIMAL
    catch err
        @warn "Mosek unavailable — skipping SDP-dependent tests" exception = err
        return false
    end
end

@testset "ppt2" begin
    include("subsystems.jl")
    include("poly2mat.jl")
    include("states.jl")
    include("ampliation.jl")
    include("gram_freedom.jl")
    include("block_positive.jl")
    include("reproduce_results.jl")
    include("scripts.jl")

    if mosek_licensed()
        include("positive_map.jl")
    else
        @testset "positive map (skipped: no Mosek)" begin
            @test_skip false
        end
    end
end
