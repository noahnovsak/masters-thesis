# PPT² final scan (seed 42, 4×4) — work, methodology, and findings

Working document for the thesis write-up. Structure and prose are complete;
`«…»` marks the exact numbers to be filled in once the run finishes (they are all
recoverable from `timings.tsv`, the per-step logs in `logs/`, and the `.jld2`
outputs / ledgers in this directory).

---

## 1. Objective

Regenerate the full 4×4 PPT² experiment under a fresh, distinct seed (42) so the
results are reproducible and independent of the earlier seed-0 data, and quantify:

1. the **detection power** of the three entanglement criteria — the level-2 DPS
   hierarchy, the PnCP **trace** witness, and the PnCP **ampliation** (positive-map)
   witness — on three families of PPT-entangled states;
2. the **PPT² conjecture** via (a) the witness SDP over PPT-map compositions for every
   generated witness, and (b) the three criteria applied to 10 000 compositions
   (ordered pairs) of states from each state family.

All runs use dimension n = m = 4, DPS level 2, tolerance 1e-8, base seed 42.

## 2. Pipeline

Run by `run_all.sh` (sequential, each step a fresh Julia process at `-t 40`),
logging wall time to `timings.tsv`. Dependency order:

| # | step | script | output | role |
|---|------|--------|--------|------|
| 1 | PnCP witness-map library | `gen_pncp.jl` | `pncp_4x4.jld2` | 10 000 PnCP maps (witnesses), seed 42 |
| 2 | witness-constructed states | `gen_witness_ppt.jl` | `witness_ppt_4x4.jld2` | one bound-entangled PPT state per witness (witness SDP) |
| 3 | random PPT, asymmetric | `compare_detection.jl` | `cmp_asym_4x4.jld2` | 10 000 detected entangled random PPT states + per-criterion scores |
| 4 | random PPT, symmetric | `compare_detection.jl --ppt-invariant` | `cmp_sym_4x4.jld2` | as (3), partial-transpose-invariant |
| 5 | detection power on witness states | `detection_power.jl` | `detection_power_witness_4x4.jld2` | trace/ampliation/DPS over the witness pool |
| 6 | PPT² via witness SDP | `gen_witness_ppt2.jl` | `witness_ppt2_4x4.jld2` | see-saw min tr(W·ampliation(ρ₁,ρ₂)) over all witnesses |
| 7–9 | PPT² on 100×100 = 10 000 pairs | `test_ppt2.jl --with-dps` | `ppt2_results/{witness,asym,sym}/` | three criteria on compositions of each state family |

## 3. Methods (definitions)

- **PnCP maps / witnesses.** `gen_pncp.jl` constructs biquadratic forms positive on the
  Segre variety but not sums of squares — Choi matrices of positive-but-not-completely-
  positive maps, i.e. entanglement witnesses. One trial = one construction attempt
  (rejected if the SOS/positivity certificate search exhausts its budget).
- **Trace criterion** (`detect_trace`): min over witnesses of `tr(W·ρ)`; entangled when
  `< −tol`. A separable state gives `tr(W·ρ) ≥ 0` for every witness.
- **Ampliation criterion** (`detect_ampliation`, **system = 1**): min over witnesses of
  `λ_min((Φ_W ⊗ I)(ρ))`; entangled when `< −tol`. System 1 is the leg dual to the
  block-positive witness — the one consistent with `tr(W·ρ) < 0`; the opposite leg
  (system 2) detects nothing for these witnesses.
- **DPS criterion** (`detect_dps`): level-2 DPS robustness (`Ket.entanglement_robustness`);
  entangled when robustness `> tol`.
- **Witness SDP** (`min_ppt_witness`): for fixed W, minimise `tr(W·ρ)` over the PPT cone;
  a negative optimum certifies a bound-entangled PPT state W detects.
- **Witness PPT² SDP** (`min_ppt2_witness`): for fixed W, see-saw minimise
  `tr(W·ampliation(ρ₁,ρ₂))` over pairs of PPT maps (16 restarts × 40 alternations); a
  negative optimum is a PPT² counterexample candidate witnessed by W. The whole-PPT-cone
  `min_ppt_witness` is its convex relaxation / lower bound.
- **Composition test** (`test_ppt2.jl`): for ordered pairs (i,j), form the composite
  `ampliation(ρᵢ, ρⱼ)` (Choi matrix of Φᵢ∘Φⱼ, itself PPT) and test it with all three
  criteria; a detected composite is a PPT² counterexample.

## 4. Witness-power comparison — framing  «important»

For comparing **witness power vs DPS**, the meaningful DPS check is the one on the
**witness-SDP states** (step 5, `detection_power.jl`). Those states are trace-detected
by construction (each by its generating witness), so the informative quantities are the
**DPS** and **ampliation** detection rates on them and the overlap (DPS∩witness,
witness-only, DPS-only). This is kept distinct from:
- the `compare_detection` random-pool analysis (steps 3–4) — detection efficacy on
  *random* PPT states, where the witnesses detect very little and DPS dominates; and
- the `test_ppt2` composition scan (steps 7–9) — the PPT² conjecture test.

## 5. Implementation & performance notes

The detection routines were benchmarked and chosen for speed (all behaviour-preserving):

- **`detect_trace`** = `findmin(W -> real(dot(W, ρ)), forms)`, using `tr(W·ρ) = real⟨W,ρ⟩`
  for real-symmetric W and Hermitian ρ. ≈ **2.5×** faster than the previous
  `tr.(forms .* Ref(ρ))` broadcast (which materialised 10 000 product matrices per call
  and drove GC pressure); numerically identical (≈1e-18).
- **`detect_ampliation`** keeps a Cholesky **early-exit** loop (sign-check via `potrf!`,
  exact `min_eig` only for the few witnesses that fire). ≈ **1.75×** faster than the
  genuine vectorised original (`eigvals.` of every witness's ampliation, commit
  `7da53c5`); identical detection.
- **`detect_dps`** drops a **negligible imaginary part** (≤1e-9 relative) before solving,
  so complex-stored-but-real states (e.g. the witness states, imag ≈1e-12) solve on the
  real PSD cone, not the double-dimension complex Hermitian one: ≈ **38×** faster
  («2.4» s vs «92» s), identical robustness; genuinely complex states are left untouched.

Shared `load_states` (layout-agnostic load + real-slice drop) lives in the `ppt2` module.

## 6. Reproducibility & run configuration

- **Seed 42** throughout (`gen_pncp`/`gen_ppt`/`compare_detection` gained a `--seed` flag
  wired to `generate_dataset`'s `seed0`; `gen_witness_ppt2` already had one). Distinct
  from the prior seed-0 datasets. `gen_witness_ppt` is deterministic (one SDP per witness).
- **Resumable**: batch generators skip completed batches; `test_ppt2` resumes from its CSV
  ledger; the driver marks completed steps and retries (3×) on failure.
- **Threads**: all steps `-t 40` with Mosek/BLAS at their defaults. (Single-thread Mosek
  pinning was measured ~1.7× *slower*; transient crashes were system pressure from other
  users, handled by resume/retry, not by throttling.) The level-2 DPS solve is
  memory-bandwidth bound at ≈1.2–1.3 solves/s, which sets the overall runtime.
- A single-threaded **warm-up** precedes each parallel loop to avoid a Julia codegen-lock
  livelock on the heavy DPS compile (otherwise the first `@threads` wave stalls > 30 min).

## 7. Results  «fill in from logs / outputs»

### 7.1 Libraries generated
- PnCP witnesses: «N» maps; acceptance rate «a»% («accepted/attempted»); time «hh:mm».
- Witness-constructed states: «N» of 10 000 witnesses certified a PPT entangled state
  («pct»%); optimum tr(W·ρ): min «·», median «·», max «·»; time «hh:mm».

### 7.2 Detection power — random PPT pools (steps 3–4)
For each pool (asymmetric, symmetric), over «N» states detected by ≥1 criterion:

| criterion | asymmetric | symmetric |
|---|---|---|
| DPS (level 2) | «n (pct%)» | «n (pct%)» |
| PnCP trace | «·» | «·» |
| PnCP ampliation (system 1) | «·» | «·» |
| PnCP (trace ∨ ampliation) | «·» | «·» |
| DPS only (PnCP missed) | «·» | «·» |
| PnCP only (DPS missed) | «·» | «·» |
| DPS ∧ PnCP | «·» | «·» |

Sampling acceptance rate ≈ «33.8»% (random PPT states that are DPS-entangled).
Expected qualitative result: PnCP witnesses detect a small fraction of random PPT
states; DPS dominates this regime.

### 7.3 Detection power — witness pool (step 5)  «witness-power result»
Over «N» witness-constructed states:

| criterion | count (pct) |
|---|---|
| DPS (level 2) | «·» |
| PnCP trace | «·» (≈100% by construction) |
| PnCP ampliation (system 1) | «·» |
| ANY | «·» |
| DPS only / PnCP only / both | «·» / «·» / «·» |

Interpretation: of the bound-entangled states the PnCP witnesses construct, DPS
independently certifies «pct»% — the witness-vs-DPS overlap.

### 7.4 PPT² via witness SDP (step 6)
Over «N» witnesses (16 restarts × 40 iters, seed 42):
- counterexample candidates (see-saw optimum `< −tol`): «n» / «N» («pct»%)
  — PPT² predicts every optimum ≥ 0; «expected: 0 counterexamples».
- see-saw optimum tr(W·composite): min «·», median «·», max «·» (boundary-distance
  distribution; how close the search gets to violating PPT²).
- counterexample witness indices (if any): «…».

### 7.5 PPT² on 10 000 compositions (steps 7–9)
For each state family, 100×100 = 10 000 ordered pairs, all three criteria (+DPS):

| pool | tested | detected entangled | trace | ampliation | dps |
|---|---|---|---|---|---|
| witness | «10000» | «n» | «·» | «·» | «·» |
| asymmetric | «10000» | «n» | «·» | «·» | «·» |
| symmetric | «10000» | «n» | «·» | «·» | «·» |

«Expected: 0 detected — every composition of PPT maps is PPT and (per PPT²) separable;
any detection is a counterexample and is saved as `result_<i>_<j>.jld2`.»

### 7.6 Timings
From `timings.tsv` (wall time per step, seed 42, -t 40):

| step | time |
|---|---|
| gen_pncp | «hh:mm:ss» |
| gen_witness_ppt | «·» |
| compare_detection_asym | «·» |
| compare_detection_sym | «·» |
| detection_power_witness | «·» |
| gen_witness_ppt2 | «·» |
| test_ppt2_witness / asym / sym | «·» / «·» / «·» |
| **total** | «·» |

## 8. Notes / caveats
- Steady-state `compare_detection` ≈ «23» min per 200-state batch at -t 40 (DPS-bound).
- DPS acceptance on random PPT states ≈ «33.8»% (240-sample estimate; refine from logs).
- The witness states are stored `ComplexF64` but real up to ~1e-12; all DPS solves run on
  the real cone (see §5).
