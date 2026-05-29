using LinearAlgebra   # eigmin in log_detection
using Base.Threads
using ProgressMeter
using JLD2            # jldopen in save_result
using ppt2            # test_ppt2, ampliation, load_batches
using ArgParse
using Ket             # partial_transpose in log_detection
using Dates
using Logging

function _parse_args()
    s = ArgParseSettings(description = "Test PPT2 over compositions of pre-generated PPT states")
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
            help = "Path to pre-generated PPT states (default: ppt_entangled_NxM.jld2)"
            arg_type = String
            default = ""
        "--forms", "-f"
            help = "Path to pre-generated PnCP forms (default: pncp_NxM.jld2)"
            arg_type = String
            default = ""
        "--output-dir", "-o"
            help = "Directory to save results and logs (default: current dir)"
            arg_type = String
            default = "."
    end
    return parse_args(s)
end

function log_detection(logger, lock, i, j, n, m, composite, robustness, min_dot, min_amp)
    eigmin_rho = eigmin(composite)
    eigmin_pt = eigmin(partial_transpose(composite, 2, [n, m]))
    Base.lock(lock) do
        with_logger(logger) do
            @info "Entanglement detected" i j eigmin_rho eigmin_pt robustness min_dot min_amp
        end
    end
end

function save_result(dir, i, j, composite, wit, dot_idx, dot_mat, amp_idx, amp_mat)
    path = joinpath(dir, "result_$(i)_$(j).jld2")
    jldopen(path, "w") do file
        file["i"] = i
        file["j"] = j
        file["state"] = composite
        file["witness"] = wit
        file["dot_idx"] = dot_idx
        file["dot_mat"] = dot_mat
        file["amp_idx"] = amp_idx
        file["amp_mat"] = amp_mat
    end
end

function run_pairs(n, m, forms, states, tol, output_dir, logger)
    N = length(states)
    log_lock = ReentrantLock()

    # all ordered pairs (i, j), including self-pairs: composition is not commutative
    @showprogress @threads for p in 0:(N * N - 1)
        i = p ÷ N + 1
        j = p % N + 1

        # compose states[i] with states[j] and run every criterion
        r = test_ppt2(states[i], states[j]; n = n, m = m, forms = forms, tol = tol, mode = :parallel)

        if r.detected
            composite = ampliation(states[i], states[j], n, m)   # only needed to log/save
            log_detection(logger, log_lock, i, j, n, m, composite, r.dps.value, r.trace.value, r.ampliation.value)
            save_result(output_dir, i, j, composite, r.dps.witness,
                        r.trace.idx, forms[r.trace.idx], r.ampliation.idx, forms[r.ampliation.idx])
        end
    end
end

function main()
    args = _parse_args()
    n = args["dim_A"]
    m = args["dim_B"]
    tol = args["tol"]
    output_dir = args["output-dir"]
    states_path = isempty(args["states"]) ? "ppt_entangled_$(n)x$(m).jld2" : args["states"]
    forms_path = isempty(args["forms"]) ? "pncp_$(n)x$(m).jld2" : args["forms"]

    mkpath(output_dir)

    log_path = joinpath(output_dir, "run_$(Dates.format(now(), "yyyymmdd_HHMMSS")).log")
    logger = SimpleLogger(open(log_path, "w"))

    println("Logging to $(log_path)")
    println("Loading states from $(states_path) ...")
    states = load_batches(states_path)
    println("Loading forms from $(forms_path) ...")
    forms = load_batches(forms_path)

    n_pairs = length(states)^2
    with_logger(logger) do
        @info "Run started" n m tol states_path forms_path output_dir n_states = length(states) n_forms = length(forms) n_pairs
    end
    println("Loaded $(length(states)) states and $(length(forms)) forms. " *
            "Testing $(n_pairs) ordered pairs on $(nthreads()) threads...")

    run_pairs(n, m, forms, states, tol, output_dir, logger)

    with_logger(logger) do
        @info "Run finished"
    end
    flush(logger.stream)
end

main()
