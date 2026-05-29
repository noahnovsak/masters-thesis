using ppt2          # rand_ppt, test_ppt2, generate_dataset, load_batches
using ArgParse

# Per accepted state we keep the state plus the raw score of every criterion;
# the detection booleans are derived from these against `tol` at summary time,
# so the threshold can be revisited without regenerating.
const DetectionResult = @NamedTuple{
    state::Matrix{Float64},
    robustness::Float64,  # DPS level-`level` robustness; entangled when > tol
    min_dot::Float64,     # min over forms of tr(form * state); entangled when < -tol
    min_amp::Float64,     # min over forms of λ_min((I⊗form)(state)); entangled when < -tol
    dot_idx::Int,         # form achieving min_dot
    amp_idx::Int,         # form achieving min_amp
}

function _parse_args()
    s = ArgParseSettings(description = "Sample random PPT states and record which criteria (PnCP witnesses vs. DPS) detect their entanglement")
    @add_arg_table! s begin
        "--total", "-t"
            help = "Total number of detected entangled states to keep"
            arg_type = Int
            default = 1000
        "--batch", "-b"
            help = "Number of kept states per batch"
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
        "--level", "-l"
            help = "DPS hierarchy level"
            arg_type = Int
            default = 2
        "--tol"
            help = "Detection tolerance"
            arg_type = Float64
            default = 1e-8
        "--ppt-invariant"
            help = "Symmetrise off-diagonal blocks so each sampled state is invariant under partial transpose"
            action = :store_true
        "--forms", "-f"
            help = "Path to pre-generated PnCP forms (default: pncp_NxM.jld2)"
            arg_type = String
            default = ""
        "--output", "-o"
            help = "Output filename (default: detection_NxM.jld2)"
            arg_type = String
            default = ""
    end
    return parse_args(s)
end

function summarize(filename, tol, level)
    results = load_batches(filename)
    N = length(results)
    N == 0 && return

    dps = [r.robustness > tol for r in results]
    trace = [r.min_dot < -tol for r in results]
    ampl = [r.min_amp < -tol for r in results]
    pncp = trace .| ampl

    pct(c) = round(100 * c / N, digits = 2)
    line(label, mask) = println("  $(rpad(label, 28))$(rpad(count(mask), 7))($(pct(count(mask)))%)")

    println("\nDetection efficacy over $(N) states (each detected by ≥1 criterion):")
    line("DPS (level $(level))", dps)
    line("PnCP trace witness", trace)
    line("PnCP ampliation", ampl)
    line("PnCP (trace OR ampliation)", pncp)
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
    ppt_invariant = args["ppt-invariant"]
    forms_path = isempty(args["forms"]) ? "pncp_$(n)x$(m).jld2" : args["forms"]
    filename = isempty(args["output"]) ? "detection_$(n)x$(m).jld2" : args["output"]

    isfile(forms_path) || error("PnCP forms not found at $(forms_path); generate them with gen_pncp.jl")
    println("Loading PnCP forms from $(forms_path) ...")
    forms = load_batches(forms_path)
    println("Loaded $(length(forms)) forms.")

    # one trial = one random PPT state, run through all criteria. The state is
    # kept when at least one criterion certifies entanglement; every score is
    # recorded so the methods can be compared afterwards.
    function trial(rng)
        state = rand_ppt(n, m; rng = rng, ppt_invariant = ppt_invariant)
        # criteria run on the state itself (not a composite), every score recorded
        r = test_ppt2(state; n = n, m = m, compose = false, forms = forms, level = level, tol = tol, mode = :parallel)

        r.detected || return nothing
        return (
            state = Matrix{Float64}(state),
            robustness = Float64(r.dps.value),
            min_dot = Float64(r.trace.value),
            min_amp = Float64(r.ampliation.value),
            dot_idx = Int(r.trace.idx),
            amp_idx = Int(r.ampliation.idx),
        )
    end

    generate_dataset(
        filename, args["total"], args["batch"], trial;
        T = DetectionResult,
        meta = Dict("dim_A" => n, "dim_B" => m, "level" => level, "tol" => tol, "ppt_invariant" => ppt_invariant, "forms" => forms_path),
        label = "detected entangled PPT states",
    )

    summarize(filename, tol, level)
end

main()
