using JLD2
using Base.Threads
using ProgressMeter
using ppt2

const TOTAL_MATRICES = 10000
const BATCH_SIZE = 1000
const N = 3
const M = 3
const FILENAME = "pncp_forms_$(N)x$(M).jld2"

function generate()
    n_batches = TOTAL_MATRICES ÷ BATCH_SIZE

    println("Starting generation on $(nthreads()) threads...")

    for b in 1:n_batches
        batch_results = Vector{Matrix{Float64}}(undef, BATCH_SIZE)

        @showprogress @threads for i in 1:BATCH_SIZE
            batch_results[i] = pncp_form(N, M)
        end

        jldopen(FILENAME, "a+") do file
            file["batch_$b"] = batch_results
        end

        println("Saved batch $b/$n_batches ($(b * BATCH_SIZE) total matrices)")
    end
end

generate()
