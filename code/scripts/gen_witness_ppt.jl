using Statistics     # median in the summary
using Base.Threads
using ProgressMeter
using JLD2            # jldopen in save
using ppt2           # min_ppt_witness, load_batches
using ArgParse

# One kept result per witness that certifies a bound entangled state: the index
# of the witness in the loaded form library, the optimal `tr(W·ρ)` (its detection
# strength, negative when entangled), and the certifying PPT state ρ itself.
const WitnessState = @NamedTuple{
    witness_idx::Int,
    value::Float64,            # min tr(W·ρ) over PPT states; entangled when < -tol
    state::Matrix{ComplexF64}, # the witnessed PPT (bound entangled) state ρ
}

function _parse_args()
    s = ArgParseSettings(description = "For each pre-generated PnCP witness W, minimise tr(W·ρ) over PPT states ρ; keep the states it certifies entangled")
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
            help = "Detection tolerance: keep witnesses with optimum below -tol"
            arg_type = Float64
            default = 1e-8
        "--limit", "-L"
            help = "Process only the first L witnesses (0 = all)"
            arg_type = Int
            default = 0
        "--forms", "-f"
            help = "Path to pre-generated PnCP forms (default: pncp_NxM.jld2)"
            arg_type = String
            default = ""
        "--output", "-o"
            help = "Output filename (default: witness_ppt_NxM.jld2)"
            arg_type = String
            default = ""
    end
    return parse_args(s)
end

function summarize(kept, N, tol)
    println("\nWitness-restricted PPT minimisation over $(N) witnesses (tol $(tol)):")
    n_det = length(kept)
    pct = N == 0 ? 0.0 : round(100 * n_det / N, digits = 2)
    println("  certified a PPT entangled state: $(n_det)/$(N) ($(pct)%)")
    if n_det > 0
        vals = [k.value for k in kept]
        println("  optimum tr(W·ρ):  min $(round(minimum(vals), sigdigits = 4)), " *
                "median $(round(median(vals), sigdigits = 4)), " *
                "max $(round(maximum(vals), sigdigits = 4))")
    end
end

function main()
    args = _parse_args()
    n = args["dim_A"]
    m = args["dim_B"]
    tol = args["tol"]
    limit = args["limit"]
    forms_path = isempty(args["forms"]) ? "pncp_$(n)x$(m).jld2" : args["forms"]
    filename = isempty(args["output"]) ? "witness_ppt_$(n)x$(m).jld2" : args["output"]

    isfile(forms_path) || error("PnCP forms not found at $(forms_path); generate them with gen_pncp.jl")
    println("Loading PnCP forms from $(forms_path) ...")
    forms = load_batches(forms_path)
    N = limit > 0 ? min(limit, length(forms)) : length(forms)
    println("Loaded $(length(forms)) forms; minimising over PPT states for $(N) on $(nthreads()) threads...")

    # one trial = one witness: minimise tr(W·ρ) over the PPT cone. The witness
    # certifies a bound entangled state exactly when the optimum drops below -tol.
    results = Vector{Union{Nothing,WitnessState}}(undef, N)
    @showprogress @threads for i in 1:N
        r = min_ppt_witness(forms[i], n, m; tol = tol)
        results[i] = r.detected ?
            (witness_idx = i, value = Float64(r.value), state = Matrix{ComplexF64}(r.state)) :
            nothing
    end

    kept = WitnessState[r for r in results if r !== nothing]
    jldopen(filename, "w") do file
        file["batch_1"] = kept
        file["meta/dim_A"] = n
        file["meta/dim_B"] = m
        file["meta/tol"] = tol
        file["meta/forms"] = forms_path
        file["meta/n_witnesses"] = N
    end
    println("Saved $(length(kept)) witnessed states to $(filename).")

    summarize(kept, N, tol)
end

main()
