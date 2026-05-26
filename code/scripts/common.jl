using JLD2
using Base.Threads
using Random
using ProgressMeter

# Shared orchestration for the generation/test scripts: resumable, reproducible,
# multithreaded batch generation with per-batch trial/success statistics.

batch_id_of(key) = parse(Int, split(key, "_")[2])

"""
    completed_batches(filename) -> Set{Int}

Ids of the `batch_<id>` groups already present in `filename` (empty if missing).
"""
function completed_batches(filename)
    isfile(filename) || return Set{Int}()
    jldopen(filename, "r") do file
        Set(batch_id_of(k) for k in keys(file) if startswith(k, "batch_"))
    end
end

"""
    load_batches(path) -> Vector{Matrix{Float64}}

Concatenate every `batch_<id>` group (in id order), skipping `meta/*` entries.
"""
function load_batches(path)
    jldopen(path, "r") do file
        ks = sort([k for k in keys(file) if startswith(k, "batch_")]; by = batch_id_of)
        isempty(ks) ? Matrix{Float64}[] : reduce(vcat, file[k] for k in ks)
    end
end

"""
    batch_counts(filename, batch_ids) -> (attempted, accepted)

Sum the stored trial/success counters for the given batches.
"""
function batch_counts(filename, batch_ids)
    attempted = 0
    accepted = 0
    (!isfile(filename) || isempty(batch_ids)) && return attempted, accepted
    jldopen(filename, "r") do file
        for id in batch_ids
            ak = "meta/batch_$(id)_attempted"
            ck = "meta/batch_$(id)_accepted"
            attempted += haskey(file, ak) ? file[ak] : 0
            accepted += haskey(file, ck) ? file[ck] : 0
        end
    end
    return attempted, accepted
end

function write_meta!(filename, meta)
    jldopen(filename, "a+") do file
        for (k, v) in meta
            mk = "meta/$(k)"
            haskey(file, mk) || (file[mk] = v)
        end
    end
end

"""
    sample_batch(trial, target, seed_base; T) -> (values, attempted)

Collect `target` accepted results. Candidate `c = 1, 2, …` is evaluated with
`trial(Xoshiro(seed_base + c))`, which returns the value to keep on success or
`nothing` on rejection. Candidates are evaluated in parallel waves; the kept
results are the `target` successes with the smallest candidate indices, and
`attempted` is the index of the last kept success.

Because each candidate's outcome is fixed by its seed and the kept set is the
lowest-index successes, the result is identical regardless of thread count or
wave sizing — the run is reproducible and thread-safe.
"""
function sample_batch(trial, target::Int, seed_base::Int; T::Type = Matrix{Float64})
    successes = Dict{Int,T}()
    next_idx = 1
    wave = max(target, nthreads())

    while length(successes) < target
        base = next_idx
        local_res = Vector{Union{Nothing,T}}(undef, wave)
        @threads for t in 1:wave
            local_res[t] = trial(Xoshiro(seed_base + base + t - 1))
        end
        for t in 1:wave
            r = local_res[t]
            r === nothing || (successes[base + t - 1] = r)
        end
        next_idx += wave

        if length(successes) < target
            evaluated = next_idx - 1
            rate = max(length(successes), 1) / evaluated
            remaining = target - length(successes)
            wave = max(nthreads(), ceil(Int, remaining / rate * 1.2))
        end
    end

    kept = sort!(collect(keys(successes)))[1:target]
    return T[successes[i] for i in kept], kept[end]
end

"""
    generate_dataset(filename, total, batch_size, trial; kwargs...)

Generate `total` accepted results in batches of `batch_size`, resuming from any
batches already in `filename`. `trial(rng)::Union{Nothing,T}` produces one
candidate. Each batch stores its results plus `meta/batch_<id>_attempted` and
`meta/batch_<id>_accepted`; entries in `meta` are written once.

Keyword arguments:
  - `seed0`   base RNG seed offset (default 0)
  - `stride`  seed span reserved per batch; must exceed the trials any batch
              needs (default 10^9)
  - `T`       element type produced by `trial` (default `Matrix{Float64}`)
  - `meta`    dataset-level metadata to record once
  - `label`   noun used in progress output (default "results")
"""
function generate_dataset(
    filename, total, batch_size, trial;
    seed0::Int = 0, stride::Int = 10^9, T::Type = Matrix{Float64},
    meta = Dict(), label::String = "results",
)
    total % batch_size == 0 || error("total ($total) must be divisible by batch ($batch_size)")
    n_batches = total ÷ batch_size
    done = completed_batches(filename)

    done_attempted, done_accepted = batch_counts(filename, done)
    if !isempty(done)
        rate = done_attempted > 0 ? 100 * done_accepted / done_attempted : 0.0
        println("Resuming: $(length(done))/$(n_batches) batches complete, " *
                "$(done_accepted)/$(done_attempted) accepted ($(round(rate, digits = 2))%)")
    end

    write_meta!(filename, meta)

    println("Generating on $(nthreads()) threads: target $(total) $(label), " *
            "$(batch_size) per batch, output $(filename)")

    total_attempted = done_attempted
    total_accepted = done_accepted
    progress = Progress(n_batches; desc = "Batches: ")
    update!(progress, length(done))

    for batch_id in 1:n_batches
        if batch_id in done
            continue
        end
        seed_base = seed0 + (batch_id - 1) * stride
        values, attempted = sample_batch(trial, batch_size, seed_base; T = T)

        jldopen(filename, "a+") do file
            file["batch_$(batch_id)"] = values
            file["meta/batch_$(batch_id)_attempted"] = attempted
            file["meta/batch_$(batch_id)_accepted"] = batch_size
        end

        total_attempted += attempted
        total_accepted += batch_size
        brate = 100 * batch_size / attempted
        trate = total_attempted > 0 ? 100 * total_accepted / total_attempted : 0.0
        next!(progress)
        println("Saved batch $(batch_id)/$(n_batches): $(batch_size)/$(attempted) accepted " *
                "($(round(brate, digits = 2))%). Running total: $(total_accepted)/$(total_attempted) " *
                "($(round(trate, digits = 2))%).")
    end

    frate = total_attempted > 0 ? 100 * total_accepted / total_attempted : 0.0
    println("Done. Accepted $(total_accepted)/$(total_attempted) ($(round(frate, digits = 2))%).")
end
