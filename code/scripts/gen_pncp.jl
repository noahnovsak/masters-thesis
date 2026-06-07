using ppt2          # pncp_mat, generate_dataset
using ArgParse

function _parse_args()
    s = ArgParseSettings(description = "Generate positive-but-not-completely-positive (PnCP) maps as matrices")
    @add_arg_table! s begin
        "--total", "-t"
            help = "Total number of maps to generate"
            arg_type = Int
            default = 1000
        "--batch", "-b"
            help = "Number of maps per batch"
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
        "--seed"
            help = "Base RNG seed (passed to generate_dataset as seed0; changes the dataset)"
            arg_type = Int
            default = 0
        "--output", "-o"
            help = "Output filename (default: pncp_NxM.jld2)"
            arg_type = String
            default = ""
    end
    return parse_args(s)
end

function main()
    args = _parse_args()
    n = args["dim_A"]
    m = args["dim_B"]
    seed = args["seed"]
    filename = isempty(args["output"]) ? "pncp_$(n)x$(m).jld2" : args["output"]

    # one trial = one construction attempt; it fails (nothing) when the
    # retry budget in pncp_mat is exhausted without a positive certificate.
    trial(rng) = pncp_mat(n, m; rng = rng)

    generate_dataset(
        filename, args["total"], args["batch"], trial;
        seed0 = seed,
        meta = Dict("dim_A" => n, "dim_B" => m, "seed" => seed),
        label = "PnCP maps",
    )
end

main()
