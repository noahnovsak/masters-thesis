# TRACE cross-detection of the witness-generated PPT entangled states.
#
# witness_ppt_4x4.jld2 holds, for every PnCP witness W_i, the bound entangled
# state ρ_i the SDP built so that tr(W_i·ρ_i) < 0 (witness i detects state i by
# the trace test, by construction). Here we test every state against the WHOLE
# witness library: T[f,s] = real tr(W_f · ρ_s), via one matrix product, and ask
# whether witnesses detect states OTHER than the one they generated.

using ppt2, JLD2, LinearAlgebra, Base.Threads, Statistics, Printf

const DATADIR = "/Users/noah/dev/masters-thesis/code/data"
const TOL = 1e-8

forms = load_batches(joinpath(DATADIR, "pncp_4x4.jld2"))
ws    = load_batches(joinpath(DATADIR, "witness_ppt_4x4.jld2"))
F, S = length(forms), length(ws)
gen_idx = Int[w.witness_idx for w in ws]
gen_val = Float64[w.value for w in ws]
@printf("Loaded %d witnesses, %d states.\n", F, S)

# T[f,s] = real tr(W_f·ρ_s) = Σ W_f[i,j]·Re ρ_s[i,j]  (W_f real-sym, ρ_s Herm.)
FM = Matrix{Float64}(undef, F, 256)
@threads for f in 1:F; @views FM[f, :] .= vec(forms[f]); end
SM = Matrix{Float64}(undef, 256, S)
@threads for s in 1:S; @views SM[:, s] .= vec(real(ws[s].state)); end
T = FM * SM

trace_self   = Float64[T[gen_idx[s], s] for s in 1:S]
trace_min    = Vector{Float64}(undef, S)
trace_argmin = Vector{Int}(undef, S)
trace_ndet   = Vector{Int}(undef, S)
@threads for s in 1:S
    col = @view T[:, s]
    v, i = findmin(col)
    trace_min[s] = v; trace_argmin[s] = i
    trace_ndet[s] = count(<(-TOL), col)
end
wtr_ndet = Vector{Int}(undef, F)
@threads for f in 1:F; wtr_ndet[f] = count(<(-TOL), @view T[f, :]); end

hr() = println("─"^64)
println(); hr()
@printf("TRACE CROSS-DETECTION  %d states × %d witnesses  (tol=%.0e)\n", S, F, TOL)
hr()
@printf("Sanity: detected by own witness: %d/%d   max|T[gen,s]−SDPval|=%.1e\n",
        count(trace_self .< -TOL), S, maximum(abs.(trace_self .- gen_val)))
tot = sum(wtr_ndet)
@printf("\nTotal (witness,state) detections: %d of %d\n", tot, S*F)
@printf("  = %d self + %d FOREIGN   (avg %.2f detecting witnesses per state)\n",
        S, tot - S, tot / S)
@printf("States detected by ≥1 FOREIGN witness:        %d/%d (%.1f%%)\n",
        count(trace_ndet .> 1), S, 100*count(trace_ndet .> 1)/S)
@printf("States whose STRONGEST witness ≠ its generator: %d/%d (%.1f%%)\n",
        count(trace_argmin .!= gen_idx), S, 100*count(trace_argmin .!= gen_idx)/S)
@printf("Witnesses detecting ≥1 FOREIGN state:          %d/%d (%.1f%%)\n",
        count(wtr_ndet .> 1), F, 100*count(wtr_ndet .> 1)/F)
@printf("  states/witness: min %d, median %.1f, mean %.2f, max %d\n",
        minimum(wtr_ndet), median(wtr_ndet), mean(wtr_ndet), maximum(wtr_ndet))
@printf("  detecting-witnesses/state: min %d, median %.1f, mean %.2f, max %d\n",
        minimum(trace_ndet), median(trace_ndet), mean(trace_ndet), maximum(trace_ndet))
bf = argmax(wtr_ndet); bs = argmax(trace_ndet)
@printf("Most prolific witness W_%d detects %d states; most-detected state ρ_%d hit by %d witnesses\n",
        bf, wtr_ndet[bf], bs, trace_ndet[bs])
hr()

jldopen(joinpath(DATADIR, "cross_trace_4x4.jld2"), "w") do f
    f["gen_idx"]=gen_idx; f["gen_val"]=gen_val
    f["trace_self"]=trace_self; f["trace_min"]=trace_min
    f["trace_argmin"]=trace_argmin; f["trace_ndet"]=trace_ndet; f["wtr_ndet"]=wtr_ndet
    f["meta/tol"]=TOL
end
println("Saved cross_trace_4x4.jld2")
