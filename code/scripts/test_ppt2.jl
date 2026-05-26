using LinearAlgebra
using JLD2
using Base.Threads
using Random
using ProgressMeter
using MosekTools
using ppt2
using ArgParse
using Ket
using Dates
using Logging

function _parse_args()
    s = ArgParseSettings(description="Test PPT2 entanglement detection over random PPT states")
    @add_arg_table! s begin
        "--dim_A", "-n"
            help = "Dimension of subspace A"
            arg_type = Int
            default = 4
        "--dim_B", "-m"
            help = "Dimension of subspace B"
            arg_type = Int
            default = 4
        "--trials", "-t"
            help = "Number of trials"
            arg_type = Int
            default = 10000
        "--tol"
            help = "Tolerance for entanglement detection"
            arg_type = Float64
            default = 1e-8
        "--forms", "-f"
            help = "Path to precomputed PNCP forms file"
            arg_type = String
            default = ""
        "--output-dir", "-o"
            help = "Directory to save results and logs (default: current dir)"
            arg_type = String
            default = "."
    end
    return parse_args(s)
end

function load_forms(path::String)
    jldopen(path, "r") do file
        vcat([file[k] for k in keys(file)]...)
    end
end

function log_detection(logger, i, state, robustness, min_dot, min_amp)
    eigmin_rho = eigmin(state)
    eigmin_pt  = eigmin(partial_transpose(state, 2))
    with_logger(logger) do
        @info "Entanglement detected in trial $i" eigmin_rho eigmin_pt robustness min_dot min_amp
    end
end

function save_result(dir, i, state, wit, dot_idx, dot_mat, amp_idx, amp_mat)
    path = joinpath(dir, "result_$(i).jld2")
    jldopen(path, "w") do file
        file["state"]   = state
        file["witness"] = wit
        file["dot_idx"] = dot_idx
        file["dot_mat"] = dot_mat
        file["amp_idx"] = amp_idx
        file["amp_mat"] = amp_mat
    end
end

function test_ppt2(n, m, forms, n_trials, tol, output_dir, logger)
    @showprogress @threads for i in 1:n_trials
        rng_i     = Xoshiro(i)
        ppt       = rand_ppt(n, m; rng=rng_i)
        composite = ampliation(ppt, ppt, n, m)

        trc, trc_i = findmin(tr.(forms .* Ref(composite)))
        amp, amp_i = findmin(minimum.(real.(eigvals.(ampliation.(forms, Ref(composite), n, m)))))
        ro, wit    = entanglement_robustness(composite, [n, m], 2; solver=Mosek.Optimizer)

        if ro > tol || trc < -tol || amp < -tol
            log_detection(logger, i, composite, ro, trc, amp)
            save_result(output_dir, i, composite, wit, trc_i, forms[trc_i], amp_i, forms[amp_i])
        end
    end
end

function main()
    args       = _parse_args()
    n          = args["dim_A"]
    m          = args["dim_B"]
    n_trials   = args["trials"]
    tol        = args["tol"]
    output_dir = args["output-dir"]
    forms_path = isempty(args["forms"]) ? "pncp_forms_$(n)x$(m).jld2" : args["forms"]

    mkpath(output_dir)

    log_path = joinpath(output_dir, "run_$(Dates.format(now(), "yyyymmdd_HHMMSS")).log")
    logger   = SimpleLogger(open(log_path, "w"))

    with_logger(logger) do
        @info "Run started" n m n_trials tol forms_path output_dir
    end

    println("Logging to $log_path")
    println("Loading forms from $forms_path ...")
    forms = load_forms(forms_path)
    println("Loaded $(length(forms)) forms. Running $n_trials trials on $(nthreads()) threads...")

    test_ppt2(n, m, forms, n_trials, tol, output_dir, logger)

    with_logger(logger) do
        @info "Run finished"
    end
end

main()
