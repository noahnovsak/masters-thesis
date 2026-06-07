using LinearAlgebra   # eigmin diagnostics in the ledger
using Base.Threads
using ProgressMeter
using JLD2            # jldopen in save_result
using ppt2           # test_ppt2, ampliation, load_batches
using ArgParse
using Ket             # partial_transpose
using Printf          # ledger formatting

# The search space (every ordered pair of states) is far too large to exhaust in
# one go, so the run is incremental and resumable. Every tested composition — not
# just the detected ones — is appended to a CSV ledger; a rerun reads it back,
# skips the pairs already there, and continues. The ledger doubles as the record
# of which compositions were tried and whether entanglement could be verified.
const LEDGER_HEADER =
    "i,j,detected,trace_value,trace_idx,amp_value,amp_idx,dps_value,eigmin_rho,eigmin_pt"

function _parse_args()
    s = ArgParseSettings(description = "Test PPT2 over compositions of pre-generated PPT states, recording every attempt in a resumable ledger")
    @add_arg_table! s begin
        "--dim_A", "-n"
            help = "Dimension of subspace A"
            arg_type = Int
            default = 4
        "--dim_B", "-m"
            help = "Dimension of subspace B"
            arg_type = Int
            default = 4
        "--tol"
            help = "Tolerance for entanglement detection"
            arg_type = Float64
            default = 1e-8
        "--states", "-s"
            help = "Path to pre-generated PPT states (default: witness_ppt_NxM.jld2)"
            arg_type = String
            default = ""
        "--forms", "-f"
            help = "Path to pre-generated PnCP forms (default: pncp_NxM.jld2)"
            arg_type = String
            default = ""
        "--output-dir", "-o"
            help = "Directory for the ledger and per-detection result files"
            arg_type = String
            default = "."
        "--ledger"
            help = "Ledger filename, relative to the output dir (default: tested_<states>.csv)"
            arg_type = String
            default = ""
        "--limit", "-L"
            help = "Test at most this many new compositions this run (0 = all remaining)"
            arg_type = Int
            default = 0
        "--max-states", "-k"
            help = "Use only the first K states, capping the pool to K*K pairs (0 = all)"
            arg_type = Int
            default = 0
        "--with-dps"
            help = "Also run the level-2 DPS robustness SDP on every pair (an SDP each — much heavier than the trace/ampliation witness criteria, but tractable on the real-valued states; off by default)"
            action = :store_true
    end
    return parse_args(s)
end

# `load_states` (layout-agnostic load + real-slice drop) lives in the ppt2 module.

# ── Resumable ledger ──────────────────────────────────────────────────────────

"Set of (i, j) pairs already recorded in `path`, plus how many were detected."
function load_ledger(path)
    done = Set{Tuple{Int,Int}}()
    n_detected = 0
    isfile(path) || return done, n_detected
    open(path) do io
        readline(io)   # header
        for line in eachline(io)
            isempty(strip(line)) && continue
            f = split(line, ',')
            push!(done, (parse(Int, f[1]), parse(Int, f[2])))
            parse(Bool, f[3]) && (n_detected += 1)
        end
    end
    return done, n_detected
end

"Open `path` for appending, writing the header first if the file is new."
function open_ledger(path)
    fresh = !isfile(path)
    io = open(path, "a")
    fresh && (println(io, LEDGER_HEADER); flush(io))
    return io
end

function write_row(io, i, j, r, eigmin_rho, eigmin_pt)
    dps_value = hasproperty(r, :dps) ? r.dps.value : NaN   # NaN when --with-dps is off
    @printf(io, "%d,%d,%s,%.17g,%d,%.17g,%d,%.17g,%.17g,%.17g\n",
            i, j, r.detected,
            r.trace.value, r.trace.idx,
            r.ampliation.value, r.ampliation.idx,
            dps_value, eigmin_rho, eigmin_pt)
end

# Rich artifact kept for each detected pair: the composite, a witness, and the two
# forms that scored best under the trace / ampliation criteria. Keys match those
# read by notebooks/read_results.ipynb.
function save_result(dir, i, j, composite, r, forms)
    # `witness`: the DPS witness when available, otherwise the trace form that
    # detected the pair (itself a valid entanglement witness).
    wit = hasproperty(r, :dps) ? r.dps.witness : forms[r.trace.idx]
    jldopen(joinpath(dir, "result_$(i)_$(j).jld2"), "w") do file
        file["i"] = i
        file["j"] = j
        file["state"] = composite
        file["witness"] = wit
        file["dot_idx"] = r.trace.idx
        file["dot_mat"] = forms[r.trace.idx]
        file["amp_idx"] = r.ampliation.idx
        file["amp_mat"] = forms[r.ampliation.idx]
    end
end

function run_pairs(n, m, forms, states, tol, output_dir, ledger_path, limit, criteria)
    K = length(states)
    total = K * K
    done, prev_detected = load_ledger(ledger_path)
    remaining = total - length(done)
    target = limit > 0 ? min(limit, remaining) : remaining
    println("Ledger $(ledger_path): $(length(done))/$(total) pairs tested " *
            "($(prev_detected) detected); $(remaining) remaining.")
    if target == 0
        println("Nothing to do.")
        return (0, 0, length(done), prev_detected)
    end

    # Warm up the per-pair path single-threaded (ampliation + criteria, incl. the heavy
    # level-2 DPS compile under --with-dps, + eigmin/partial_transpose) so the first
    # @threads wave doesn't livelock on Julia's codegen lock (observed >30 min stall).
    println("Warming up the per-pair path (single-threaded compile)...")
    let comp = Matrix(ampliation(states[1], states[1], n, m))
        test_ppt2(comp; n = n, m = m, compose = false, forms = forms,
                  criteria = criteria, tol = tol, mode = :parallel)
        eigmin(Hermitian(comp))
        eigmin(Hermitian(partial_transpose(comp, 2, [n, m])))
    end

    io = open_ledger(ledger_path)
    chunk = max(nthreads(), 64)   # parallel wave; ledger is flushed after each
    processed = 0
    detected_run = 0
    p = 0                         # linear cursor over ordered pairs, i = p÷K+1, j = p%K+1
    prog = Progress(target; desc = "Pairs: ")

    try
        while p < total && processed < target
            # gather the next wave of untested pairs (skipping any already in the ledger)
            idxs = Int[]
            while length(idxs) < chunk && p < total && processed + length(idxs) < target
                i = p ÷ K + 1
                j = p % K + 1
                (i, j) in done || push!(idxs, p)
                p += 1
            end
            isempty(idxs) && break

            results = Vector{Any}(undef, length(idxs))
            @threads for t in eachindex(idxs)
                i = idxs[t] ÷ K + 1
                j = idxs[t] % K + 1
                comp = Matrix(ampliation(states[i], states[j], n, m))
                r = test_ppt2(comp; n = n, m = m, compose = false,
                              forms = forms, criteria = criteria, tol = tol, mode = :parallel)
                eigmin_rho = eigmin(Hermitian(comp))
                eigmin_pt = eigmin(Hermitian(partial_transpose(comp, 2, [n, m])))
                results[t] = (; i, j, comp, r, eigmin_rho, eigmin_pt)
            end

            # persist single-threaded: append a ledger row for every pair and a
            # rich result file for the detected ones
            for res in results
                write_row(io, res.i, res.j, res.r, res.eigmin_rho, res.eigmin_pt)
                if res.r.detected
                    detected_run += 1
                    save_result(output_dir, res.i, res.j, res.comp, res.r, forms)
                end
                push!(done, (res.i, res.j))
            end
            flush(io)
            processed += length(idxs)
            update!(prog, processed)
        end
    finally
        finish!(prog)
        close(io)
    end
    return (processed, detected_run, length(done), prev_detected + detected_run)
end

function main()
    args = _parse_args()
    n = args["dim_A"]
    m = args["dim_B"]
    tol = args["tol"]
    output_dir = args["output-dir"]
    limit = args["limit"]
    max_states = args["max-states"]
    # The two PnCP-witness criteria are cheap matrix ops; the DPS SDP is much heavier
    # (an SDP per pair), so it is opt-in (--with-dps) rather than run on every pair.
    criteria = args["with-dps"] ? (:trace, :ampliation, :dps) : (:trace, :ampliation)
    states_path = isempty(args["states"]) ? "witness_ppt_$(n)x$(m).jld2" : args["states"]
    forms_path = isempty(args["forms"]) ? "pncp_$(n)x$(m).jld2" : args["forms"]
    ledger_name = isempty(args["ledger"]) ?
        "tested_$(splitext(basename(states_path))[1]).csv" : args["ledger"]
    ledger_path = joinpath(output_dir, ledger_name)

    mkpath(output_dir)

    println("Loading states from $(states_path) ...")
    states = load_states(states_path)
    if max_states > 0 && max_states < length(states)
        states = states[1:max_states]
    end
    println("Loading forms from $(forms_path) ...")
    forms = load_batches(forms_path)

    K = length(states)
    println("Loaded $(K) states and $(length(forms)) forms. " *
            "Pool of $(K * K) ordered pairs; criteria $(join(criteria, '+')); " *
            "testing on $(nthreads()) threads (limit $(limit == 0 ? "none" : limit)).")

    processed, detected_run, total_done, total_detected =
        run_pairs(n, m, forms, states, tol, output_dir, ledger_path, limit, criteria)

    println("\nThis run: tested $(processed) compositions, $(detected_run) detected entangled.")
    println("Cumulative: $(total_done) tested, $(total_detected) detected " *
            "(of $(K * K) possible pairs).")
end

main()
