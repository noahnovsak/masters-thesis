using JLD2
using Base.Threads
using Random
using ProgressMeter
using ppt2
using ArgParse
using Ket
using MosekTools

function sample_random_ppt(n, m; rng=Random.GLOBAL_RNG)
    if isdefined(ppt2, :ran_ppt)
        return getfield(ppt2, :ran_ppt)(n, m; rng=rng)
    end
    return rand_ppt(n, m; rng=rng)
end

function _parse_args()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--total", "-t"
            help = "Total number of random PPT states to sample"
            arg_type = Int
            default = 1000
        "--batch", "-b"
            help = "Number of sampled states per batch"
            arg_type = Int
            default = 200
        "--dim_A", "-n"
            help = "Dimension of subspace A"
            arg_type = Int
            default = 4
        "--dim_B", "-m"
            help = "Dimension of subspace B"
            arg_type = Int
            default = 4
        "--tol"
            help = "Entanglement tolerance: keep states with robustness > tol"
            arg_type = Float64
            default = 1e-8
        "--output", "-o"
            help = "Output filename (default: ppt_entangled_NxM.jld2)"
            arg_type = String
            default = ""
    end
    return parse_args(s)
end

function completed_batches(filename)
    isfile(filename) || return Set{Int}()
    jldopen(filename, "r") do file
        Set(
            parse(Int, split(key, "_")[2]) for key in keys(file)
            if startswith(key, "batch_")
        )
    end
end

function batch_counts(filename, batch_ids)
    attempted = 0
    accepted = 0

    if !isfile(filename) || isempty(batch_ids)
        return attempted, accepted
    end

    jldopen(filename, "r") do file
        for batch_id in batch_ids
            attempted += haskey(file, "meta/batch_$(batch_id)_attempted") ? file["meta/batch_$(batch_id)_attempted"] : 0
            accepted += haskey(file, "meta/batch_$(batch_id)_accepted") ? file["meta/batch_$(batch_id)_accepted"] : 0
        end
    end

    return attempted, accepted
end

function ensure_metadata(filename, n, m, tol)
    jldopen(filename, "a+") do file
        if !haskey(file, "meta/tol")
            file["meta/tol"] = tol
        end
        if !haskey(file, "meta/dim_A")
            file["meta/dim_A"] = n
        end
        if !haskey(file, "meta/dim_B")
            file["meta/dim_B"] = m
        end
    end
end

function generate(total_states, batch_size, n, m, tol, filename)
    n_batches = total_states ÷ batch_size
    done = completed_batches(filename)

    done_attempted, done_accepted = batch_counts(filename, done)
    if !isempty(done)
        done_rate = done_attempted > 0 ? 100 * done_accepted / done_attempted : 0.0
        println("Resuming: $(length(done))/$(n_batches) batches already complete")
        println("Completed stats so far: $(done_accepted)/$(done_attempted) entangled ($(round(done_rate, digits=2))%)")
    end

    println("Starting generation on $(nthreads()) threads...")
    println("Config: $(total_states) accepted states target, $(batch_size) per batch, $(n)x$(m), tol=$(tol), output: $(filename)")

    total_attempted = done_attempted
    total_accepted = done_accepted

    ensure_metadata(filename, n, m, tol)

    for batch_id in 1:n_batches
        if batch_id in done
            continue
        end

        entangled_states = Vector{Matrix{Float64}}(undef, batch_size)
        batch_attempted = 0
        batch_accepted = 0

        p = ProgressUnknown(desc="Batch $(batch_id)/$(n_batches)")
        while batch_accepted < batch_size
            batch_attempted += 1
            rng_i = Xoshiro((batch_id - 1) * 10^7 + batch_attempted)
            state = sample_random_ppt(n, m; rng=rng_i)
            robustness, _ = entanglement_robustness(state, [n, m], 2; solver=Mosek.Optimizer)

            if robustness > tol
                batch_accepted += 1
                entangled_states[batch_accepted] = Matrix{Float64}(state)
                next!(p; showvalues=[(:accepted, "$(batch_accepted)/$(batch_size)"), (:attempted, batch_attempted)])
            elseif batch_attempted % 25 == 0
                next!(p; showvalues=[(:accepted, "$(batch_accepted)/$(batch_size)"), (:attempted, batch_attempted)])
            end
        end
        finish!(p)

        batch_rate = batch_attempted > 0 ? 100 * batch_accepted / batch_attempted : 0.0

        jldopen(filename, "a+") do file
            file["batch_$(batch_id)"] = entangled_states
            file["meta/batch_$(batch_id)_attempted"] = batch_attempted
            file["meta/batch_$(batch_id)_accepted"] = batch_accepted
        end

        total_attempted += batch_attempted
        total_accepted += batch_accepted
        total_rate = total_attempted > 0 ? 100 * total_accepted / total_attempted : 0.0

        println(
            "Saved batch $(batch_id)/$(n_batches): " *
            "$(batch_accepted)/$(batch_attempted) entangled ($(round(batch_rate, digits=2))%). " *
            "Running total: $(total_accepted)/$(total_attempted) ($(round(total_rate, digits=2))%)."
        )
    end

    final_rate = total_attempted > 0 ? 100 * total_accepted / total_attempted : 0.0
    println("Done. Success rate: $(total_accepted)/$(total_attempted) = $(round(final_rate, digits=2))%")
end

function main()
    args = _parse_args()

    total_states = args["total"]
    batch_size = args["batch"]
    n = args["dim_A"]
    m = args["dim_B"]
    tol = args["tol"]
    filename = isempty(args["output"]) ? "ppt_entangled_$(n)x$(m).jld2" : args["output"]

    if total_states <= 0
        error("--total must be positive")
    end
    if batch_size <= 0
        error("--batch must be positive")
    end
    if total_states % batch_size != 0
        error("--total must be divisible by --batch")
    end

    generate(total_states, batch_size, n, m, tol, filename)
end

main()
