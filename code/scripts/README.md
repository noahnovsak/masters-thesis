# Scripts

Pipeline for the PPT2 search: generate a pool of PnCP witness maps and a pool of
entangled PPT states, then test every composition of those states for
entanglement.

All scripts share [`common.jl`](common.jl), which provides resumable,
reproducible, multithreaded batch generation with per-batch trial/success
statistics.

## Setup

Run everything from the `code/` directory (the Julia project root). Install
dependencies once:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Generation and testing are multithreaded — set the thread count with `-t`:

```sh
julia --project=. -t auto scripts/<script>.jl [options]
```

The SDP steps use Mosek, so a valid Mosek license must be available.

## Pipeline

Run the generators first, then the test:

```sh
# 1. PnCP witness maps  ->  pncp_4x4.jld2
julia --project=. -t auto scripts/gen_pncp.jl --total 10000 --batch 200 -n 4 -m 4

# 2a. Entangled PPT states by random sampling  ->  ppt_entangled_4x4.jld2
julia --project=. -t auto scripts/gen_ppt.jl --total 1000 --batch 200 -n 4 -m 4 --tol 1e-8

# 2b. ...or one bound entangled state per witness  ->  witness_ppt_4x4.jld2
julia --project=. -t auto scripts/gen_witness_ppt.jl -n 4 -m 4 --tol 1e-8

# 3. Test PPT2 over ordered pairs of states, using the forms as witnesses.
#    Incremental and resumable: run a slice at a time with --limit.
julia --project=. -t auto scripts/test_ppt2.jl -n 4 -m 4 --tol 1e-8 --limit 1000
```

### `gen_pncp.jl`

Generates positive-but-not-completely-positive maps (as symmetric matrices).
One trial is one construction attempt; it is rejected when `pncp_mat` exhausts
its retry budget without a positive, non-SOS certificate.

| option | default | meaning |
| --- | --- | --- |
| `--total`, `-t` | 1000 | total maps to keep |
| `--batch`, `-b` | 200 | maps per batch |
| `--dim_A`, `-n` | 4 | dimension of subspace A |
| `--dim_B`, `-m` | 4 | dimension of subspace B |
| `--seed` | 0 | base RNG seed (passed to `generate_dataset` as `seed0`; changes the dataset) |
| `--output`, `-o` | `pncp_NxM.jld2` | output file |

### `gen_ppt.jl`

Samples random PPT states and keeps the entangled ones. One trial is one random
PPT state; it is accepted when the level-2 DPS relaxation certifies
entanglement (robustness `> --tol`).

| option | default | meaning |
| --- | --- | --- |
| `--total`, `-t` | 1000 | entangled states to keep |
| `--batch`, `-b` | 200 | states per batch |
| `--dim_A`, `-n` | 4 | dimension of subspace A |
| `--dim_B`, `-m` | 4 | dimension of subspace B |
| `--tol` | 1e-8 | keep states with robustness above this |
| `--ppt-invariant` | off | symmetrise off-diagonal blocks so each state is invariant under partial transpose |
| `--seed` | 0 | base RNG seed (passed to `generate_dataset` as `seed0`; changes the dataset) |
| `--output`, `-o` | `ppt_entangled_NxM.jld2` | output file |

### `test_ppt2.jl`

Loads the pre-generated states and forms, then for **ordered pairs** `(i, j)`
(including self-pairs — composition is not commutative) forms the composite
`ampliation(states[i], states[j])` and checks it for entanglement via three
tests: the minimum `tr(form · composite)` over all forms and the minimum
eigenvalue of `ampliation(form, composite)` over all forms — both cheap
matrix scans over the form library — plus, only with `--with-dps`, the level-2
DPS robustness SDP. A pair is detected when any active test exceeds `--tol`.

The DPS SDP is much heavier than the witness criteria — a semidefinite program
per pair — so it is **off by default**; leave it off for bulk scanning and turn it
on to cross-check specific pairs. The states are stored complex but are
real-symmetric up to SDP noise, so they are loaded as real, which keeps the DPS
solve on the real PSD cone instead of the (double-dimension) complex Hermitian one
— ≈20× faster and ≈9× less memory (~2 min vs ~47 min, ~31 GB vs ~283 GB for 64
pairs on 8 threads).

The state pool is large (`witness_ppt_4x4.jld2` alone is 10000 states, i.e. 10⁸
ordered pairs), so the run is **incremental and resumable**. Every tested pair —
detected or not — is appended as one row to a CSV ledger
`tested_<states>.csv`, recording the outcome and each criterion's score (see
[Output format](#output-format)). A rerun reads the ledger back, skips the pairs
already in it, and continues; `--limit` caps how many new pairs a single run
tests, and `--max-states` caps the pool to the first `K` states (`K*K` pairs).
Detected pairs additionally get a `result_<i>_<j>.jld2` with the composite,
witness, and best detecting forms.

Accepts either state-file layout: the bare matrices from `gen_ppt.jl` or the
`(witness_idx, value, state)` tuples from `gen_witness_ppt.jl` (the default).

| option | default | meaning |
| --- | --- | --- |
| `--dim_A`, `-n` | 4 | dimension of subspace A |
| `--dim_B`, `-m` | 4 | dimension of subspace B |
| `--tol` | 1e-8 | detection tolerance |
| `--states`, `-s` | `witness_ppt_NxM.jld2` | pre-generated PPT states |
| `--forms`, `-f` | `pncp_NxM.jld2` | pre-generated PnCP forms |
| `--output-dir`, `-o` | `.` | directory for the ledger and result files |
| `--ledger` | `tested_<states>.csv` | ledger filename (relative to the output dir) |
| `--limit`, `-L` | 0 | test at most this many new pairs this run (0 = all remaining) |
| `--max-states`, `-k` | 0 | use only the first K states (0 = all) |
| `--with-dps` | off | also run the level-2 DPS SDP per pair (an SDP each — much heavier than the witness criteria) |

The `dps_value` ledger column is `nan` for pairs tested without `--with-dps`.

### `compare_detection.jl`

Compares the efficacy of the two entanglement-detection methods — PnCP witness
maps vs. the DPS hierarchy — on random PPT states. Like `gen_ppt.jl` it samples
random PPT states, but it runs **every** criterion on each state and keeps any
state that at least one criterion flags, recording each criterion's raw score:
the minimum `tr(form · state)` and minimum eigenvalue of `(I⊗form)(state)` over
the form library (the two PnCP criteria), and the level-`--level` DPS robustness.
After generating it prints a per-criterion breakdown (including DPS-only,
PnCP-only, and overlap counts) so the methods can be compared directly. Requires
a pre-generated PnCP form library (`gen_pncp.jl`).

```sh
julia --project=. -t auto scripts/compare_detection.jl --total 1000 --batch 200 -n 4 -m 4
```

| option | default | meaning |
| --- | --- | --- |
| `--total`, `-t` | 1000 | detected entangled states to keep |
| `--batch`, `-b` | 200 | states per batch |
| `--dim_A`, `-n` | 4 | dimension of subspace A |
| `--dim_B`, `-m` | 4 | dimension of subspace B |
| `--level`, `-l` | 2 | DPS hierarchy level |
| `--tol` | 1e-8 | detection tolerance |
| `--ppt-invariant` | off | symmetrise off-diagonal blocks so each sampled state is invariant under partial transpose |
| `--seed` | 0 | base RNG seed (passed to `generate_dataset` as `seed0`; changes the dataset) |
| `--forms`, `-f` | `pncp_NxM.jld2` | pre-generated PnCP forms |
| `--output`, `-o` | `detection_NxM.jld2` | output file |

### `detection_power.jl`

Measures the detection power of the three criteria — DPS, the PnCP trace witness, and
the PnCP ampliation (`system=1`) — on a **pre-generated** state library, rather than
sampling its own like `compare_detection.jl`. Loads an existing state file (any layout,
via the shared `load_states`) and runs every criterion on every state, then prints the
same per-criterion breakdown (DPS-only / PnCP-only / overlap). Use it to score a fixed
pool such as the witness-constructed states from `gen_witness_ppt.jl`. Requires a
pre-generated PnCP form library.

```sh
julia --project=. -t auto scripts/detection_power.jl -n 4 -m 4 --level 2 \
    -s witness_ppt_4x4.jld2 -f pncp_4x4.jld2
```

| option | default | meaning |
| --- | --- | --- |
| `--dim_A`, `-n` | 4 | dimension of subspace A |
| `--dim_B`, `-m` | 4 | dimension of subspace B |
| `--level`, `-l` | 2 | DPS hierarchy level |
| `--tol` | 1e-8 | detection tolerance |
| `--states`, `-s` | — | pre-generated state library to scan (required) |
| `--forms`, `-f` | `pncp_NxM.jld2` | pre-generated PnCP forms |
| `--output`, `-o` | `detection_power_NxM.jld2` | per-state scores output |

### `gen_witness_ppt.jl`

A third entanglement test, dual to the trace witness. Instead of fixing a state
and scanning the form library for the most negative `tr(W · ρ)`, it fixes each
witness `W` and solves the SDP

```
minimise   tr(W · ρ)
subject to ρ ⪰ 0,  tr(ρ) = 1,  ρ^Γ ⪰ 0   (Hermitian, PSD, unit-trace, PPT)
```

over the whole PPT cone (`ppt2.min_ppt_witness`). Since every separable state
gives `tr(W · ρ) ≥ 0`, a negative optimum certifies the minimiser as a PPT
**entangled** (bound entangled) state that `W` detects — so the run both measures
each witness's detection strength and emits the bound entangled states it
witnesses. Iterates over a pre-generated form library in parallel; requires
`gen_pncp.jl` output. The companion notebook
[`notebooks/sdp_witness_ppt.ipynb`](../notebooks/sdp_witness_ppt.ipynb) drives the
same SDP interactively.

```sh
julia --project=. -t auto scripts/gen_witness_ppt.jl -n 4 -m 4 --tol 1e-8
```

| option | default | meaning |
| --- | --- | --- |
| `--dim_A`, `-n` | 4 | dimension of subspace A |
| `--dim_B`, `-m` | 4 | dimension of subspace B |
| `--tol` | 1e-8 | keep witnesses whose optimum is below `-tol` |
| `--limit`, `-L` | 0 | process only the first L witnesses (0 = all) |
| `--forms`, `-f` | `pncp_NxM.jld2` | pre-generated PnCP forms |
| `--output`, `-o` | `witness_ppt_NxM.jld2` | output file |

### `gen_witness_ppt2.jl`

Sharpens `gen_witness_ppt.jl` from the whole PPT cone down to the **composition
manifold** — the PPT² setting itself. Instead of minimising `tr(W · ρ)` over an
arbitrary PPT state `ρ`, it minimises `tr(W · ampliation(ρ1, ρ2))` over *pairs* of
PPT maps `ρ1, ρ2` (`ppt2.min_ppt2_witness`):

```
minimise   tr(W · composite),   composite = ampliation(ρ1, ρ2)
subject to ρ1, ρ2 ⪰ 0,  tr = 1,  ρ^Γ ⪰ 0   (each factor Hermitian, PSD, unit-trace, PPT)
```

`ampliation` is bilinear in `(ρ1, ρ2)`, so this is not a single SDP; it is solved
by **see-saw** (freeze one factor, optimise the other, alternate) from
`--restarts` random PPT starts. A composition of PPT maps is itself PPT, so a
negative optimum exhibits a PPT *and* entangled composition — a **PPT²
counterexample** witnessed by `W`. Because see-saw only finds a local optimum, a
non-negative result is supporting evidence, not proof; `gen_witness_ppt.jl` (the
whole PPT cone) is the convex relaxation / lower bound. Each witness `i` is seeded
`Xoshiro(--seed + i)`, so the run is reproducible and independent of thread count.
Requires `gen_pncp.jl` output; the companion notebook
[`notebooks/sdp_witness_ppt.ipynb`](../notebooks/sdp_witness_ppt.ipynb) drives the
underlying SDPs interactively.

```sh
julia --project=. -t auto scripts/gen_witness_ppt2.jl -n 4 -m 4 --tol 1e-8 --restarts 16 --max_iter 40
```

| option | default | meaning |
| --- | --- | --- |
| `--dim_A`, `-n` | 4 | dimension of subspace A |
| `--dim_B`, `-m` | 4 | dimension of subspace B |
| `--tol` | 1e-8 | keep witnesses whose optimum is below `-tol` |
| `--restarts` | 16 | see-saw random restarts per witness |
| `--max_iter` | 40 | max alternating SDP steps per restart |
| `--seed` | 0 | base RNG seed (witness `i` uses `Xoshiro(seed + i)`) |
| `--limit`, `-L` | 0 | process only the first L witnesses (0 = all) |
| `--forms`, `-f` | `pncp_NxM.jld2` | pre-generated PnCP forms |
| `--output`, `-o` | `witness_ppt2_NxM.jld2` | output file |

### `cross_trace.jl` and `cross_ampl.jl`

The witness-construction SDP builds, for each PnCP witness `W_i`, exactly the one
bound entangled state `ρ_i` that `W_i` detects by the trace test. These two scripts
ask how far each witness reaches *beyond its own* state: they cross-evaluate the
whole witness library against the whole witness-derived state pool (10⁴ × 10⁴ ≈ 10⁸
pairs) and report how many foreign states each witness detects.

- **`cross_trace.jl`** computes the trace functional `T[f,s] = Re tr(W_f · ρ_s)` for
  every (witness, state) pair as a single matrix product and counts detections
  (`T < -tol`) per witness and per state.
- **`cross_ampl.jl`** runs the stronger, nonlinear ampliation test
  `λ_min((Φ_{W_f} ⊗ I)(ρ_s))` over the same grid, with an allocation-free,
  block-checkpointed inner loop that resumes after a kill.

Unlike the other scripts these take **no command-line options**; they read the
pre-generated `pncp_4x4.jld2` and `witness_ppt_4x4.jld2` from a data directory and
write their results (`cross_trace_4x4.jld2`, `cross_ampl_*_4x4.jld2`) back to it.
The directory defaults to `results/` and is overridable with the `DATADIR`
environment variable:

```sh
DATADIR=results julia --project=. -t auto scripts/cross_trace.jl
DATADIR=results julia --project=. -t auto scripts/cross_ampl.jl
```

The finding (thesis §Results) is that each witness detects essentially only the
state extracted from it — the library behaves as a collection of single-state
detectors.

## Output format

Generated `.jld2` files store data under `batch_<id>` keys and statistics under
`meta/*` (`dim_A`, `dim_B`, `tol`, and per-batch `batch_<id>_attempted` /
`batch_<id>_accepted` counts). The generators store a `Vector{Matrix}` per
batch; `compare_detection.jl` stores a `Vector` of named tuples
`(state, robustness, min_dot, min_amp, dot_idx, amp_idx)` instead, and
`gen_witness_ppt.jl` a `Vector` of `(witness_idx, value, state)` (the
witness's index in the form library, its optimum `tr(W · ρ)`, and the certified
PPT entangled state as a `Matrix{ComplexF64}`). `gen_witness_ppt2.jl` stores a
`Vector` of `(witness_idx, value, rho1, rho2)` for the detected compositions only
(the two PPT factors whose composition is the counterexample candidate), plus
`meta/all_values` — every witness's see-saw optimum, for the boundary-distance
distribution — alongside the `meta/restarts`, `meta/max_iter`, and `meta/seed`
settings.

`test_ppt2.jl` instead writes a plain-CSV **ledger** `tested_<states>.csv` — one
row per tested ordered pair, header
`i,j,detected,trace_value,trace_idx,amp_value,amp_idx,dps_value,eigmin_rho,eigmin_pt`
— which is the resume state (pairs already present are skipped) and the full
record of which compositions were tried and whether entanglement was verified.
Each detected pair also gets a `result_<i>_<j>.jld2` holding the composite
`state`, the DPS `witness`, and the best trace/ampliation forms with their
indices (`dot_idx`/`dot_mat`, `amp_idx`/`amp_mat`).

## Resumability and reproducibility

- **Resumable:** rerun a generator with the same `--output` to fill in only the
  missing batches; completed batches are detected and skipped.
- **Reproducible:** every candidate is seeded deterministically, and a batch
  keeps the lowest-index successes, so a given configuration produces the same
  dataset regardless of thread count. Changing `--batch` changes the seed layout
  and therefore the dataset.
