using JLD2
using Base.Threads
using Random
using ProgressMeter
using ppt2

const TOTAL_MATRICES = 1000
const BATCH_SIZE = 200
const N = 4
const M = 4
const RNG = Xoshiro(0)
const FILENAME = "pncp_forms_$(N)x$(M).jld2"

function generate()
    n_batches = TOTAL_MATRICES ÷ BATCH_SIZE

    println("Starting generation on $(nthreads()) threads...")

    for b in 1:n_batches
        batch_results = Vector{Matrix{Float64}}(undef, BATCH_SIZE)

        @showprogress @threads for i in 1:BATCH_SIZE
            batch_results[i] = pncp_mat(N, M, RNG)
        end

        jldopen(FILENAME, "a+") do file
            file["batch_$b"] = batch_results
        end

        println("Saved batch $b/$n_batches ($(b * BATCH_SIZE) total matrices)")
    end
end

generate()
