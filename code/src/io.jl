# All dataset I/O: the scripts generate batched `.jld2` files, the notebooks
# read them back. Both sides share these.

using JLD2: jldopen
using Base.Threads: @threads, nthreads
using ProgressMeter: Progress, update!, next!

# ── Readers ──────────────────────────────────────────────────────────────────

batch_id_of(key) = parse(Int, split(key, "_")[2])

"""
    load_batches(path) -> Vector

Concatenate every `batch_<id>` group in `path`, in id order (`meta/*` skipped).
This is the ordering `dot_idx`/`amp_idx` and a composite's `i`/`j` index into.
"""
function load_batches(path)
    jldopen(path, "r") do file
        ks = sort([k for k in keys(file) if startswith(k, "batch_")]; by = batch_id_of)
        isempty(ks) ? [] : reduce(vcat, file[k] for k in ks)
    end
end

"""
    load_meta(path) -> Dict{String,Any}

The `meta` group as a Dict (dataset values plus per-batch attempted/accepted
counters); empty when the file carries no metadata.
"""
function load_meta(path)
    jldopen(path, "r") do file
        haskey(file, "meta") || return Dict{String,Any}()
        g = file["meta"]
        Dict{String,Any}(k => g[k] for k in keys(g))
    end
end

# ── Writers: resumable, reproducible, multithreaded batch generation ──────────

"Ids of the `batch_<id>` groups already in `filename` (empty if missing)."
function completed_batches(filename)
    isfile(filename) || return Set{Int}()
    jldopen(filename, "r") do file
        Set(batch_id_of(k) for k in keys(file) if startswith(k, "batch_"))
    end
end

"Sum the stored attempted/accepted counters for the given batches."
function batch_counts(filename, batch_ids)
    attempted = accepted = 0
    (!isfile(filename) || isempty(batch_ids)) && return attempted, accepted
    jldopen(filename, "r") do file
        for id in batch_ids
            ak, ck = "meta/batch_$(id)_attempted", "meta/batch_$(id)_accepted"
            attempted += haskey(file, ak) ? file[ak] : 0
            accepted += haskey(file, ck) ? file[ck] : 0
        end
    end
    return attempted, accepted
end

"Write each `meta` entry once; existing keys are left untouched."
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

Collect `target` accepted results. Candidate `c` is `trial(Xoshiro(seed_base + c))`,
returning the value to keep or `nothing`. Candidates run in parallel waves; the
kept set is the `target` lowest-index successes, so the result is independent of
thread count. `attempted` is the last kept success's index.
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
            rate = max(length(successes), 1) / (next_idx - 1)
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
batches already in `filename`. `trial(rng)::Union{Nothing,T}` is one candidate.
Each batch stores its results plus `meta/batch_<id>_attempted`/`_accepted`.

Keywords: `seed0` (base seed offset), `stride` (seed span per batch, must exceed
any batch's trials), `T` (element type), `meta` (dataset metadata, written once),
`label` (noun for progress output).
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
        batch_id in done && continue
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
