# AMPLIATION cross-detection of the witness-generated PPT entangled states.
#
# For every witness W_f and state ρ_s:  A[f,s] = λ_min( (I_n ⊗ Φ_{W_f})(ρ_s) ).
# Entanglement is certified when A[f,s] < -tol. This nonlinear positive-map test
# is stronger than the trace functional, so we ask whether it detects FOREIGN
# states (f ≠ gen_idx[s]) that the trace test (which found none) misses.
#
# Allocation-free inner loop (per-thread buffers + mul! + in-place Hermitian eig)
# and block checkpointing to cross_ampl_4x4.jld2 so the run RESUMES after a kill.

using ppt2, JLD2, LinearAlgebra, Base.Threads, Statistics, Printf

BLAS.set_num_threads(1)
const DATADIR = "/Users/noah/dev/masters-thesis/code/data"
# Map acts on subsystem A (= the subsystem dual to the trace witness). We feed the
# SWAPPED state SW·ρ·SW into the standard (I⊗Φ) ampliation, since
# (I⊗Φ)(SW·ρ·SW) = (Φ⊗I)(ρ). The plain (I⊗Φ)(ρ) version detects nothing (wrong
# subsystem) — see /tmp/cross_ampl.log and the verify scripts.
const SUB1 = true
const CKPT = joinpath(DATADIR, SUB1 ? "cross_ampl_sub1_4x4.jld2" : "cross_ampl_4x4.jld2")
const TOL = 1e-8
const n, m, r = 4, 4, 4
const BLK = 500                                    # states per checkpoint block

forms = load_batches(joinpath(DATADIR, "pncp_4x4.jld2"))
ws    = load_batches(joinpath(DATADIR, "witness_ppt_4x4.jld2"))
F, S = length(forms), length(ws)
gen_idx = Int[w.witness_idx for w in ws]
gen_val = Float64[w.value for w in ws]

# natural rep of each W_f (as ComplexF64 so mul! hits zgemm) and regrouped ρ_s
natrep(W) = Matrix{ComplexF64}(reshape(PermutedDimsArray(reshape(W, r,m,r,m), (1,3,2,4)), r*r, m*m))
bmat(ρ)   = Matrix{ComplexF64}(reshape(PermutedDimsArray(reshape(ρ, m,n,m,n), (1,3,2,4)), m*m, n*n))
natAs = Vector{Matrix{ComplexF64}}(undef, F); @threads for f in 1:F; natAs[f]=natrep(forms[f]); end
const SW = SUB1 ? Matrix{ComplexF64}(ppt2.swap(n)) : Matrix{ComplexF64}(I, n*m, n*m)
Bms   = Vector{Matrix{ComplexF64}}(undef, S)
@threads for s in 1:S; Bms[s]=bmat(SW * Matrix{ComplexF64}(ws[s].state) * SW); end

# Build M = (I⊗Φ)(ρ) into preallocated Mbuf from natA (form) and Bm (state).
@inline function buildM!(Cbuf, Mbuf, natA, Bm)
    mul!(Cbuf, natA, Bm)                           # C[(k,l),(α,β)]
    @inbounds for β in 1:n, l in 1:r, α in 1:n, k in 1:r
        Mbuf[k+(α-1)*r, l+(β-1)*r] = Cbuf[k+(l-1)*r, α+(β-1)*n]
    end
    return Mbuf
end

# Fast path: is λ_min(M) ≥ −TOL?  Test M+TOL·I ≻ 0 by in-place Cholesky (potrf!,
# ~6× cheaper than an eigendecomposition and allocation-free). Destroys Mbuf.
@inline function amp_notdetected!(Cbuf, Mbuf, natA, Bm)
    buildM!(Cbuf, Mbuf, natA, Bm)
    @inbounds for i in 1:n*r; Mbuf[i, i] += TOL; end
    _, info = LinearAlgebra.LAPACK.potrf!('U', Mbuf)
    return info == 0                               # true ⟺ not detected
end

# Exact smallest eigenvalue (only used for the rare detected pairs + self pair).
@inline amp_exact!(Cbuf, Mbuf, natA, Bm) = eigvals!(Hermitian(buildM!(Cbuf, Mbuf, natA, Bm)))[1]

# ── resume from checkpoint ────────────────────────────────────────────────────
amp_self  = fill(NaN, S); amp_minv = fill(NaN, S)
amp_argmin= zeros(Int, S); amp_ndet = fill(-1, S)
wamp_ndet = zeros(Int, F); done_upto = 0
if isfile(CKPT)
    jldopen(CKPT, "r") do f
        done_upto = f["done_upto"]
        amp_self[1:done_upto]  = f["amp_self"][1:done_upto]
        amp_minv[1:done_upto]  = f["amp_min"][1:done_upto]
        amp_argmin[1:done_upto]= f["amp_argmin"][1:done_upto]
        amp_ndet[1:done_upto]  = f["amp_ndet"][1:done_upto]
        wamp_ndet .= f["wamp_ndet"]
    end
    @printf("Resuming from checkpoint: %d/%d states already done.\n", done_upto, S)
end
@printf("Ampliation scan: %d witnesses × %d states, %d threads, from state %d.\n",
        F, S, nthreads(), done_upto+1)

function scan!(amp_self, amp_minv, amp_argmin, amp_ndet, wamp_resumed,
               Bms, natAs, gen_idx, gen_val, S, F, done_upto)
    nt = nthreads()
    wamp_local = [zeros(Int, F) for _ in 1:nt]
    t0 = time()
    for b in (done_upto ÷ BLK + 1):cld(S, BLK)
        rng = ((b-1)*BLK+1):min(b*BLK, S)
        first(rng) <= done_upto && (rng = (done_upto+1):last(rng))
        isempty(rng) && continue
        sub = collect(Iterators.partition(rng, cld(length(rng), nt)))
        @threads for ci in eachindex(sub)
            Cbuf = Matrix{ComplexF64}(undef, r*r, n*n)
            Mbuf = Matrix{ComplexF64}(undef, n*r, n*r)
            lc = wamp_local[ci]
            for s in sub[ci]
                Bm = Bms[s]; g = gen_idx[s]
                best = Inf; bestf = 0; ndet = 0
                for f in 1:F
                    amp_notdetected!(Cbuf, Mbuf, natAs[f], Bm) && continue
                    v = amp_exact!(Cbuf, Mbuf, natAs[f], Bm)      # exact λ_min, detected pair only
                    ndet += 1; lc[f] += 1
                    v < best && (best = v; bestf = f)
                end
                self = amp_exact!(Cbuf, Mbuf, natAs[g], Bm)       # self strength (one pair, cheap)
                bestf == 0 && (best = self; bestf = g)            # nothing detected → report self
                amp_self[s]=self; amp_minv[s]=best; amp_argmin[s]=bestf; amp_ndet[s]=ndet
            end
        end
        done_upto = last(rng)
        wsum = reduce(+, wamp_local) .+ wamp_resumed
        jldopen(CKPT, "w") do f
            f["done_upto"]=done_upto; f["amp_self"]=amp_self; f["amp_min"]=amp_minv
            f["amp_argmin"]=amp_argmin; f["amp_ndet"]=amp_ndet; f["wamp_ndet"]=wsum
            f["gen_idx"]=gen_idx; f["gen_val"]=gen_val; f["meta/tol"]=TOL
        end
        @printf("  %d/%d states (%.0f%%, %.1f min)  foreign-amp-dets so far: %d\n",
                done_upto, S, 100*done_upto/S, (time()-t0)/60,
                sum(amp_ndet[i] - (amp_self[i]<-TOL) for i in 1:done_upto if amp_ndet[i]>=0)); flush(stdout)
    end
    return reduce(+, wamp_local) .+ wamp_resumed
end

wamp_ndet = scan!(amp_self, amp_minv, amp_argmin, amp_ndet, wamp_ndet,
                  Bms, natAs, gen_idx, gen_val, S, F, done_upto)

# ── report ────────────────────────────────────────────────────────────────────
hr() = println("─"^64)
amp_self_det = count(<(-TOL), amp_self)
amp_any_det  = count(<(-TOL), amp_minv)
amp_cross    = count(s -> amp_ndet[s] - (amp_self[s]<-TOL) > 0, 1:S)
tot = sum(wamp_ndet)
println(); hr()
@printf("AMPLIATION CROSS-DETECTION  %d states × %d witnesses (tol=%.0e)\n", S, F, TOL)
hr()
@printf("total (f,s) detections:                 %d of %d\n", tot, S*F)
@printf("states detected by their OWN witness:   %d/%d (%.1f%%)\n", amp_self_det, S, 100*amp_self_det/S)
@printf("states detected by ANY witness:         %d/%d (%.1f%%)\n", amp_any_det, S, 100*amp_any_det/S)
@printf("states detected by ≥1 FOREIGN witness:  %d/%d (%.1f%%)\n", amp_cross, S, 100*amp_cross/S)
@printf("witnesses detecting ≥1 state (ampl.):   %d/%d (%.1f%%)\n", count(>(0),wamp_ndet), F, 100*count(>(0),wamp_ndet)/F)
@printf("states/witness (ampl.): min %d, median %.1f, mean %.2f, max %d\n",
        minimum(wamp_ndet), median(wamp_ndet), mean(wamp_ndet), maximum(wamp_ndet))
@printf("self λ_min(I⊗Φ_s)(ρ_s): min %.3e, median %.3e, max %.3e\n",
        minimum(amp_self), median(amp_self), maximum(amp_self))
hr()
println("Done. Checkpoint/results in $(CKPT)")
