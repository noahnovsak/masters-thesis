# PPT² final scan (seed 42, 4×4) — work, methodology, and findings

Working document for the thesis write-up. The run is **complete** and all numbers are
filled in (§7); the raw artifacts are in `timings.tsv`, the per-step logs in `logs/`,
the `.jld2` outputs, and the CSV ledgers in this directory. See `SUMMARY.md` for the
condensed numbers-and-timings table.

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
| 3 | random PPT, asymmetric | `compare_detection.jl` | `cmp_asym_4x4.jld2` | 5 000 detected entangled random PPT states + per-criterion scores |
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

## 4. Witness-power comparison — framing

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
  (2.5 s vs 92 s), identical robustness; genuinely complex states are left untouched.

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

## 7. Results

### 7.1 Libraries generated
- PnCP witnesses: **10 000** maps; acceptance rate **100.0%** (every construction attempt
  yielded a valid certificate); time 39 m 12 s.
- Witness-constructed states: **10 000 / 10 000** witnesses certified a PPT entangled state
  (**100.0%**); optimum tr(W·ρ): min −0.03484, median −5.292e-5, max −6.057e-8; time 3 m 34 s.

### 7.2 Detection power — random PPT pools (steps 3–4)
Over 5 000 states per pool (each detected by ≥1 criterion):

| criterion | asymmetric | symmetric |
|---|---|---|
| DPS (level 2) | 5000 (100.0%) | 5000 (100.0%) |
| PnCP trace | **0 (0.0%)** | **0 (0.0%)** |
| PnCP ampliation (system 1) | **0 (0.0%)** | **0 (0.0%)** |
| PnCP (trace ∨ ampliation) | 0 (0.0%) | 0 (0.0%) |
| DPS only (PnCP missed) | 5000 (100.0%) | 5000 (100.0%) |
| PnCP only (DPS missed) | 0 (0.0%) | 0 (0.0%) |
| DPS ∧ PnCP | 0 (0.0%) | 0 (0.0%) |

Sampling acceptance rate ≈ 33.8% (random PPT states that are DPS-entangled). **Result:
the 10 000-witness PnCP library detects none of the random PPT-entangled states; DPS
certifies all of them. The witnesses do not generalise to generic random PPT states.**

### 7.3 Detection power — witness pool (step 5) — witness-power result
Over 10 000 witness-constructed states:

| criterion | count (pct) |
|---|---|
| DPS (level 2) | 10000 (100.0%) |
| PnCP trace | 10000 (100.0%) (by construction) |
| PnCP ampliation (system 1) | 10000 (100.0%) |
| ANY | 10000 (100.0%) |
| DPS only / PnCP only / both | 0 / 0 / **10000 (100.0%)** |

**Interpretation: of the bound-entangled states the PnCP witnesses construct, DPS
independently certifies 100% — complete witness-vs-DPS overlap. Ampliation (system 1)
also fires on all 10 000, tracking the trace witness exactly. So on the witnesses' own
states the two methods are equally powerful; the contrast with §7.2 (0% on random states)
shows the witnesses are specialised to the states they construct.**

### 7.4 PPT² via witness SDP (step 6)
Over 10 000 witnesses (16 restarts × 40 iters, seed 42), time 1 h 33 m 30 s:
- counterexample candidates (see-saw optimum `< −tol`): **0 / 10 000 (0.0%)**.
- see-saw optimum tr(W·composite): min 1.557e-16, median 2.985e-10, max 2.237e-9 — every
  optimum ≥ 0, consistent with PPT² (the search gets as close as ~1e-16 to the boundary
  but never crosses it).
- counterexample witness indices: none.

### 7.5 PPT² on 10 000 compositions (steps 7–9)
For each state family, 100×100 = 10 000 ordered pairs, all three criteria (+DPS):

| pool | tested | detected entangled | trace | ampliation | dps |
|---|---|---|---|---|---|
| witness | 10000 | **0** | 0 | 0 | 0 |
| asymmetric | 10000 | **0** | 0 | 0 | 0 |
| symmetric | 10000 | **0** | 0 | 0 | 0 |

**0 detected across all 30 000 compositions — no PPT² counterexample by any criterion.**
Combined with §7.4 (0/10 000 via the witness SDP), no method found a single PPT² violation
at 4×4: strong numerical support for the conjecture.

### 7.6 Timings
From `timings.tsv` (wall time per step, seed 42, -t 40):

| step | time |
|---|---|
| gen_pncp | 0:39:12 |
| gen_witness_ppt | 0:03:34 |
| compare_detection_asym (5000) | (assembled across sessions; ≈ sym below) |
| compare_detection_sym (5000) | 4:31:04 |
| detection_power_witness | 6:04:56 |
| gen_witness_ppt2 | 1:33:30 |
| test_ppt2_witness / asym / sym | 3:32:15 / 4:36:39 / 4:15:34 |
| **total (steps with clean timing)** | **≈ 25.2 h** |

## 8. Notes / caveats
- The random-PPT pools (steps 3–4) were capped at **5 000** states each (down from
  10 000): running DPS on every sampled state is the long pole, and checking that many
  states for the detection-power comparison is unnecessary (the pair test caps the pool
  to 100, and the witness pool carries the witness-vs-DPS comparison). The asymmetric
  pool was truncated from an in-progress 7 800-state run to exactly 5 000 (25 batches),
  which equals what a fresh `--total 5000` seed-42 run produces. Scaling back up later
  is a `--total` change away (the generators resume).
- Steady-state `compare_detection` ≈ 11 min per 200-state batch at -t 40 on an idle box
  (sym: 25 batches in 4 h 31 m); ≈ 23 min/batch earlier under other-user load. DPS-bound.
- DPS acceptance on random PPT states ≈ 33.8% (240-sample estimate).
- The witness states are stored `ComplexF64` but real up to ~1e-12; all DPS solves run on
  the real cone (see §5).
