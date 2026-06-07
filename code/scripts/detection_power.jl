using ppt2          # test_ppt2, load_batches, load_states-style helpers
using ArgParse
using Base.Threads
using ProgressMeter
using JLD2

# Measure the detection power of the three entanglement criteria — DPS, the PnCP
# trace witness, and the PnCP ampliation (system=1, the leg dual to the witness) —
# on a *pre-generated* library of PPT states. Unlike compare_detection.jl, which
# samples its own states and keeps the detected ones, this loads an existing state
# file (e.g. the witness-constructed bound entangled states from gen_witness_ppt.jl)
# and runs every criterion on every state, so all three numbers are measured on the
# same fixed pool. Accepts any state-file layout (bare matrices, the
# (witness_idx,value,state) tuples from gen_witness_ppt.jl, or the named tuples from
# compare_detection.jl) via the same real-part handling as test_ppt2.jl.

const Score = @NamedTuple{
    robustness::Float64,  # DPS level-`level` robustness; entangled when > tol
    min_dot::Float64,     # min over forms of tr(form * state); entangled when < -tol
    min_amp::Float64,     # min over forms of λ_min((Φ_form ⊗ I)(state)); entangled when < -tol
    dot_idx::Int,
    amp_idx::Int,
}

function _parse_args()
    s = ArgParseSettings(description = "Measure DPS / trace / ampliation detection power on a pre-generated PPT-state library")
    @add_arg_table! s begin
        "--dim_A", "-n"
            arg_type = Int
            default = 4
        "--dim_B", "-m"
            arg_type = Int
            default = 4
        "--level", "-l"
            help = "DPS hierarchy level"
            arg_type = Int
            default = 2
        "--tol"
            help = "Detection tolerance"
            arg_type = Float64
            default = 1e-8
        "--states", "-s"
            help = "Path to the pre-generated PPT states to scan"
            arg_type = String
            default = ""
        "--forms", "-f"
            help = "Path to pre-generated PnCP forms (default: pncp_NxM.jld2)"
            arg_type = String
            default = ""
        "--output", "-o"
            help = "Output filename for the per-state scores (default: detection_power_NxM.jld2)"
            arg_type = String
            default = ""
    end
    return parse_args(s)
end

# `load_states` (layout-agnostic load + real-slice drop) lives in the ppt2 module.

function summarize(scores, tol, level, label)
    N = length(scores)
    N == 0 && return
    dps   = [s.robustness > tol for s in scores]
    trace = [s.min_dot < -tol for s in scores]
    ampl  = [s.min_amp < -tol for s in scores]
    pncp  = trace .| ampl
    any_  = dps .| pncp

    pct(c) = round(100 * c / N, digits = 2)
    line(lbl, mask) = println("  $(rpad(lbl, 28))$(rpad(count(mask), 7))($(pct(count(mask)))%)")

    println("\nDetection power over $(N) states from $(label):")
    line("DPS (level $(level))", dps)
    line("PnCP trace witness", trace)
    line("PnCP ampliation (system=1)", ampl)
    line("PnCP (trace OR ampliation)", pncp)
    line("ANY criterion", any_)
    println("  " * "-"^36)
    line("DPS only (PnCP missed)", dps .& .!pncp)
    line("PnCP only (DPS missed)", pncp .& .!dps)
    line("DPS and PnCP", dps .& pncp)
end

function main()
    args = _parse_args()
    n = args["dim_A"]
    m = args["dim_B"]
    level = args["level"]
    tol = args["tol"]
    isempty(args["states"]) && error("--states is required")
    states_path = args["states"]
    forms_path = isempty(args["forms"]) ? "pncp_$(n)x$(m).jld2" : args["forms"]
    filename = isempty(args["output"]) ? "detection_power_$(n)x$(m).jld2" : args["output"]

    isfile(states_path) || error("states not found at $(states_path)")
    isfile(forms_path) || error("PnCP forms not found at $(forms_path); generate them with gen_pncp.jl")

    println("Loading states from $(states_path) ...")
    states = load_states(states_path)
    println("Loading PnCP forms from $(forms_path) ...")
    forms = load_batches(forms_path)
    N = length(states)
    println("Scanning $(N) states against $(length(forms)) forms " *
            "(trace + ampliation(system=1) + DPS level $(level)) on $(nthreads()) threads...")

    # Warm up the criterion path (esp. the heavy level-2 DPS compile) single-threaded
    # so the first @threads wave doesn't livelock on Julia's codegen lock.
    println("Warming up the criterion path (single-threaded compile)...")
    test_ppt2(states[1]; n = n, m = m, compose = false, forms = forms,
              criteria = (:trace, :ampliation, :dps), level = level, tol = tol, mode = :parallel)

    scores = Vector{Score}(undef, N)
    @showprogress @threads for i in 1:N
        r = test_ppt2(states[i]; n = n, m = m, compose = false, forms = forms,
                      criteria = (:trace, :ampliation, :dps), level = level, tol = tol, mode = :parallel)
        scores[i] = (
            robustness = Float64(r.dps.value),
            min_dot = Float64(r.trace.value),
            min_amp = Float64(r.ampliation.value),
            dot_idx = Int(r.trace.idx),
            amp_idx = Int(r.ampliation.idx),
        )
    end

    jldopen(filename, "w") do file
        file["batch_1"] = scores
        file["meta/dim_A"] = n
        file["meta/dim_B"] = m
        file["meta/level"] = level
        file["meta/tol"] = tol
        file["meta/states"] = states_path
        file["meta/forms"] = forms_path
    end
    println("Saved $(N) per-state scores to $(filename).")

    summarize(scores, tol, level, basename(states_path))
end

main()
