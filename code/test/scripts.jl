using Test
using JLD2
using Random

# `common.jl` is the shared batch-generation engine behind every script
# (gen_pncp / gen_ppt / test_ppt2 / compare_detection).
include(joinpath(@__DIR__, "..", "scripts", "common.jl"))

# A deterministic trial: candidate `c` (seeded `Xoshiro(seed_base + c)`) is
# accepted when its first random draw is below `thresh`, returning that draw as a
# 1×1 matrix. Re-deriving the per-seed outcome lets the tests predict exactly
# which candidates the scheduler must keep.
make_trial(thresh) = rng -> (h = rand(rng); h < thresh ? fill(h, 1, 1) : nothing)
seed_value(seed) = rand(Xoshiro(seed))

@testset "batch_id_of" begin
    @test batch_id_of("batch_1") == 1
    @test batch_id_of("batch_42") == 42
end

@testset "sample_batch keeps the lowest-index successes" begin
    thresh, target, seed_base = 0.3, 5, 1000
    values, attempted = sample_batch(make_trial(thresh), target, seed_base)

    expected = Int[]                         # recompute kept set from first principles
    idx = 0
    while length(expected) < target
        idx += 1
        seed_value(seed_base + idx) < thresh && push!(expected, idx)
    end
    @test length(values) == target
    @test attempted == expected[end]         # "attempted = last kept success"
    @test [v[1, 1] for v in values] ≈ [seed_value(seed_base + i) for i in expected]
end

@testset "sample_batch is reproducible" begin
    trial = make_trial(0.4)
    a_vals, a_att = sample_batch(trial, 8, 555)
    b_vals, b_att = sample_batch(trial, 8, 555)
    @test a_att == b_att
    @test a_vals == b_vals
end

@testset "generate_dataset round-trip + resume + reproducibility" begin
    mktempdir() do dir
        trial = make_trial(0.5)
        f1 = joinpath(dir, "ds.jld2")

        generate_dataset(f1, 2, 2, trial; meta=Dict("dim_A" => 2))   # 1 batch of 2
        @test completed_batches(f1) == Set(1)
        first_batch = load_batches(f1)
        @test length(first_batch) == 2

        generate_dataset(f1, 4, 2, trial; meta=Dict("dim_A" => 2))   # resume to 2 batches
        @test completed_batches(f1) == Set([1, 2])
        all_vals = load_batches(f1)
        @test length(all_vals) == 4
        @test all_vals[1:2] == first_batch       # batch 1 preserved on resume

        attempted, accepted = batch_counts(f1, [1, 2])
        @test accepted == 4                      # accepted == batch size
        @test attempted ≥ accepted

        f2 = joinpath(dir, "ds2.jld2")           # same config ⟹ identical dataset
        generate_dataset(f2, 4, 2, trial)
        @test load_batches(f2) == all_vals
    end
end

@testset "generate_dataset requires divisible total" begin
    mktempdir() do dir
        @test_throws ErrorException generate_dataset(joinpath(dir, "x.jld2"), 5, 2, make_trial(0.5))
    end
end

@testset "write_meta! writes each key once" begin
    mktempdir() do dir
        f = joinpath(dir, "m.jld2")
        write_meta!(f, Dict("a" => 1, "b" => 2))
        write_meta!(f, Dict("a" => 99))          # must not overwrite an existing key
        jldopen(f, "r") do file
            @test file["meta/a"] == 1
            @test file["meta/b"] == 2
        end
    end
end

@testset "missing-file accessors return empties" begin
    f = joinpath(mktempdir(), "nope.jld2")
    @test completed_batches(f) == Set{Int}()
    @test batch_counts(f, [1, 2]) == (0, 0)
end
