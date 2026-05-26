using Random
using ppt2
using ArgParse
using Ket
using MosekTools

include(joinpath(@__DIR__, "common.jl"))

function _parse_args()
    s = ArgParseSettings(description = "Sample random PPT states and keep the entangled ones")
    @add_arg_table! s begin
        "--total", "-t"
            help = "Total number of entangled PPT states to keep"
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

function main()
    args = _parse_args()
    n = args["dim_A"]
    m = args["dim_B"]
    tol = args["tol"]
    filename = isempty(args["output"]) ? "ppt_entangled_$(n)x$(m).jld2" : args["output"]

    # one trial = one random PPT state; accepted only when the level-2 DPS
    # relaxation certifies entanglement (robustness above tol).
    function trial(rng)
        state = rand_ppt(n, m; rng = rng)
        robustness, _ = entanglement_robustness(state, [n, m], 2; solver = Mosek.Optimizer)
        return robustness > tol ? Matrix{Float64}(state) : nothing
    end

    generate_dataset(
        filename, args["total"], args["batch"], trial;
        meta = Dict("dim_A" => n, "dim_B" => m, "tol" => tol),
        label = "entangled PPT states",
    )
end

main()
