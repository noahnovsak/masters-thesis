using Random
using JLD2          # jldopen for manual batch writes
using ProgressMeter # Progress, update!, next!
using ppt2          # gram_freedom, load_batches, completed_batches, write_meta!
using ArgParse

# Expand a pool of PnCP forms (e.g. pncp_4x4.jld2) into equivalent Gram-matrix
# representations. The Gram matrix of a biquadratic form is unique only up to the
# space L of matrices vanishing on the real Segre variety (`gram_freedom`), so for
# each source form M0 we emit `count` alternatives M0 + Σ λ·N that represent the
# *same* polynomial but differ off the Segre variety — a family of candidate
# witnesses sharing one separable boundary.

function _parse_args()
    s = ArgParseSettings(description = "Expand each PnCP form into several equivalent (asymmetric) Gram-matrix representations of the same polynomial")
    @add_arg_table! s begin
        "--count", "-c"
            help = "Number of Gram representations to emit per source form"
            arg_type = Int
            default = 10
        "--scale"
            help = "Standard deviation of the Gram-freedom coefficients λ"
            arg_type = Float64
            default = 1.0
        "--dim_A", "-n"
            help = "Dimension of subspace A"
            arg_type = Int
            default = 4
        "--dim_B", "-m"
            help = "Dimension of subspace B"
            arg_type = Int
            default = 4
        "--seed"
            help = "Base RNG seed"
            arg_type = Int
            default = 0
        "--input", "-i"
            help = "Input PnCP forms (default: pncp_NxM.jld2)"
            arg_type = String
            default = ""
        "--output", "-o"
            help = "Output filename (default: pncp_NxM_asym.jld2)"
            arg_type = String
            default = ""
    end
    return parse_args(s)
end

function main()
    args = _parse_args()
    n = args["dim_A"]
    m = args["dim_B"]
    count = args["count"]
    scale = args["scale"]
    seed0 = args["seed"]
    input = isempty(args["input"]) ? "pncp_$(n)x$(m).jld2" : args["input"]
    output = isempty(args["output"]) ? "pncp_$(n)x$(m)_asym.jld2" : args["output"]

    isfile(input) || error("input forms not found at $(input); generate them with gen_pncp.jl")
    forms = load_batches(input)
    N = length(forms)
    N == 0 && error("no forms found in $(input)")

    # Basis of the Gram freedom L (dim = binomial(n,2)·binomial(m,2)). Adding any
    # element of span(B) to a Gram matrix leaves the polynomial unchanged.
    B = gram_freedom(n, m)

    write_meta!(output, Dict(
        "dim_A" => n, "dim_B" => m, "count" => count, "scale" => scale,
        "seed" => seed0, "source" => input, "source_forms" => N,
    ))

    done = completed_batches(output)
    println("Expanding $(N) forms from $(input) into $(count) Gram representations each " *
            "(L has $(length(B)) generators) → $(output)")
    isempty(done) || println("Resuming: $(length(done))/$(N) source forms already done")

    progress = Progress(N; desc = "Forms: ")
    update!(progress, length(done))

    # One source form → one batch of `count` representations. Each batch is seeded
    # by its source index alone (Xoshiro(seed0 + f)), so the dataset is identical
    # regardless of evaluation order and a rerun fills only the missing batches.
    for f in 1:N
        if f in done
            next!(progress)
            continue
        end

        M0 = forms[f]
        rng = Xoshiro(seed0 + f)
        reps = Matrix{Float64}[M0 + scale * sum(randn(rng) * Nα for Nα in B) for _ in 1:count]

        jldopen(output, "a+") do file
            file["batch_$(f)"] = reps
            file["meta/batch_$(f)_attempted"] = count
            file["meta/batch_$(f)_accepted"] = count
        end

        next!(progress)
    end

    println("Done. Wrote $(count) representations for $(N) forms ($(N * count) total) to $(output).")
end

main()
