using JLD2
using Base.Threads
using Random
using ProgressMeter
using ppt2
using ArgParse

function _parse_args()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--total", "-t"
            help = "Total number of matrices to generate"
            arg_type = Int
            default = 1000
        "--batch", "-b"
            help = "Number of matrices per batch"
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
        "--output", "-o"
            help = "Output filename (default: pncp_forms_NxM.jld2)"
            arg_type = String
            default = ""
    end
    return parse_args(s)
end

function completed_batches(filename)
    isfile(filename) || return Set{Int}()
    jldopen(filename, "r") do file
        Set(parse(Int, split(key, "_")[2]) for key in keys(file))
    end
end

function generate(total_matrices, batch_size, N, M, filename)
    n_batches = total_matrices ÷ batch_size
    done      = completed_batches(filename)

    if !isempty(done)
        println("Resuming: $(length(done))/$(n_batches) batches already complete, skipping $(sort(collect(done)))")
    end

    println("Starting generation on $(nthreads()) threads...")
    println("Config: $(total_matrices) matrices, $(batch_size) per batch, $(N)x$(M), output: $(filename)")

    for batch_id in 1:n_batches
        if batch_id in done
            continue
        end

        batch_results = Vector{Matrix{Float64}}(undef, batch_size)

        @showprogress @threads for i in 1:batch_size
            rng_i = Xoshiro(batch_id * batch_size + i)  # unique, stable seed per element
            batch_results[i] = pncp_mat(N, M; rng=rng_i)
        end

        jldopen(filename, "a+") do file
            file["batch_$(batch_id)"] = batch_results
        end

        println("Saved batch $batch_id/$n_batches ($(batch_id * batch_size) total matrices)")
    end
end

function main()
    args = _parse_args()

    total_matrices = args["total"]
    batch_size     = args["batch"]
    N              = args["dim_A"]
    M              = args["dim_B"]
    filename       = isempty(args["output"]) ? "pncp_forms_$(N)x$(M).jld2" : args["output"]

    generate(total_matrices, batch_size, N, M, filename)
end

main()
