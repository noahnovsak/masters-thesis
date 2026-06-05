using Statistics     # median in the summary
using Base.Threads
using Random          # per-witness Xoshiro seeding: thread-safe and reproducible
using ProgressMeter
using JLD2            # jldopen in save
using ppt2           # min_ppt2_witness, load_batches
using ArgParse

# One kept result per witness whose see-saw found a *composition* it detects as
# entangled — a candidate PPT² counterexample. Stores the witness's index in the
# form library, the optimal tr(W·composite) (negative when the composition is
# detected entangled), and the two PPT factors ρ1, ρ2 whose composition
# `ampliation(ρ1, ρ2)` attains it (the composite itself is derivable, so not saved).
const WitnessComposition = @NamedTuple{
    witness_idx::Int,
    value::Float64,             # min tr(W·composite) over compositions; counterexample when < -tol
    rho1::Matrix{ComplexF64},   # first PPT factor  Φ_1
    rho2::Matrix{ComplexF64},   # second PPT factor Φ_2
}

function _parse_args()
    s = ArgParseSettings(description = "For each pre-generated PnCP witness W, see-saw minimise tr(W·ampliation(ρ1,ρ2)) over pairs of PPT maps; keep any composition it detects as entangled (a PPT² counterexample)")
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
        "--restarts"
            help = "See-saw random restarts per witness"
            arg_type = Int
            default = 16
        "--max_iter"
            help = "Max alternating SDP steps per restart"
            arg_type = Int
            default = 40
        "--seed"
            help = "Base RNG seed (witness i uses Xoshiro(seed + i))"
            arg_type = Int
            default = 0
        "--limit", "-L"
            help = "Process only the first L witnesses (0 = all)"
            arg_type = Int
            default = 0
        "--forms", "-f"
            help = "Path to pre-generated PnCP forms (default: pncp_NxM.jld2)"
            arg_type = String
            default = ""
        "--output", "-o"
            help = "Output filename (default: witness_ppt2_NxM.jld2)"
            arg_type = String
            default = ""
    end
    return parse_args(s)
end

function summarize(kept, all_values, N, tol)
    println("\nSee-saw PPT² witness search over $(N) witnesses (tol $(tol)):")
    n_det = length(kept)
    pct = N == 0 ? 0.0 : round(100 * n_det / N, digits = 2)
    println("  detected a composition (PPT² counterexample candidate): $(n_det)/$(N) ($(pct)%)")
    if !isempty(all_values)
        println("  see-saw optimum tr(W·composite):  min $(round(minimum(all_values), sigdigits = 4)), " *
                "median $(round(median(all_values), sigdigits = 4)), " *
                "max $(round(maximum(all_values), sigdigits = 4))")
        println("  (PPT² predicts every optimum ≥ 0; a value below -tol is a counterexample candidate.)")
    end
    n_det > 0 && println("  ⚠ counterexample candidates at witness indices: $([k.witness_idx for k in kept])")
end

function main()
    args = _parse_args()
    n = args["dim_A"]
    m = args["dim_B"]
    tol = args["tol"]
    restarts = args["restarts"]
    max_iter = args["max_iter"]
    seed = args["seed"]
    limit = args["limit"]
    forms_path = isempty(args["forms"]) ? "pncp_$(n)x$(m).jld2" : args["forms"]
    filename = isempty(args["output"]) ? "witness_ppt2_$(n)x$(m).jld2" : args["output"]

    isfile(forms_path) || error("PnCP forms not found at $(forms_path); generate them with gen_pncp.jl")
    println("Loading PnCP forms from $(forms_path) ...")
    forms = load_batches(forms_path)
    N = limit > 0 ? min(limit, length(forms)) : length(forms)
    println("Loaded $(length(forms)) forms; see-saw over PPT-map compositions for $(N) witnesses " *
            "($(restarts) restarts × $(max_iter) iters) on $(nthreads()) threads...")

    # one trial = one witness: see-saw minimise tr(W·composite) over compositions
    # of two PPT maps. Witness i is seeded Xoshiro(seed + i), so the result depends
    # only on the witness and its seed — reproducible and independent of thread
    # count. A composition is detected (and kept as a PPT² counterexample candidate)
    # when the optimum drops below -tol. `all_values` records every witness's
    # optimum for the boundary-distance distribution, even the (expected) non-hits.
    all_values = Vector{Float64}(undef, N)
    results = Vector{Union{Nothing,WitnessComposition}}(undef, N)
    @showprogress @threads for i in 1:N
        r = min_ppt2_witness(forms[i], n, m; restarts = restarts, max_iter = max_iter,
                             tol = tol, rng = Xoshiro(seed + i))
        all_values[i] = Float64(r.value)
        results[i] = r.detected ?
            (witness_idx = i, value = Float64(r.value),
             rho1 = Matrix{ComplexF64}(r.ρ1), rho2 = Matrix{ComplexF64}(r.ρ2)) :
            nothing
    end

    kept = WitnessComposition[r for r in results if r !== nothing]
    jldopen(filename, "w") do file
        file["batch_1"] = kept
        file["meta/dim_A"] = n
        file["meta/dim_B"] = m
        file["meta/tol"] = tol
        file["meta/restarts"] = restarts
        file["meta/max_iter"] = max_iter
        file["meta/seed"] = seed
        file["meta/forms"] = forms_path
        file["meta/n_witnesses"] = N
        file["meta/all_values"] = all_values
    end
    println("Saved $(length(kept)) counterexample candidates to $(filename).")

    summarize(kept, all_values, N, tol)
end

main()
