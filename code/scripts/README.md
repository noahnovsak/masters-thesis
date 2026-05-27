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

Run the two generators first, then the test:

```sh
# 1. PnCP witness maps  ->  pncp_4x4.jld2
julia --project=. -t auto scripts/gen_pncp.jl --total 10000 --batch 200 -n 4 -m 4

# 2. Entangled PPT states  ->  ppt_entangled_4x4.jld2
julia --project=. -t auto scripts/gen_ppt.jl --total 1000 --batch 200 -n 4 -m 4 --tol 1e-8

# 3. Test PPT2 over all ordered pairs of states, using the forms as witnesses
julia --project=. -t auto scripts/test_ppt2.jl -n 4 -m 4 --tol 1e-8
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
| `--output`, `-o` | `ppt_entangled_NxM.jld2` | output file |

### `test_ppt2.jl`

Loads the pre-generated states and forms, then for **every ordered pair**
`(i, j)` (including self-pairs — composition is not commutative) forms the
composite `ampliation(states[i], states[j])` and checks it for entanglement via
three tests: the minimum `tr(form · composite)` over all forms, the minimum
eigenvalue of `ampliation(form, composite)` over all forms, and the level-2 DPS
robustness SDP. A pair is flagged when any test exceeds `--tol`; its composite,
witness, and best detecting forms are written to `result_<i>_<j>.jld2`.

| option | default | meaning |
| --- | --- | --- |
| `--dim_A`, `-n` | 4 | dimension of subspace A |
| `--dim_B`, `-m` | 4 | dimension of subspace B |
| `--tol` | 1e-8 | detection tolerance |
| `--states`, `-s` | `ppt_entangled_NxM.jld2` | pre-generated PPT states |
| `--forms`, `-f` | `pncp_NxM.jld2` | pre-generated PnCP forms |
| `--output-dir`, `-o` | `.` | directory for results and the run log |

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
| `--forms`, `-f` | `pncp_NxM.jld2` | pre-generated PnCP forms |
| `--output`, `-o` | `detection_NxM.jld2` | output file |

## Output format

Generated `.jld2` files store data under `batch_<id>` keys and statistics under
`meta/*` (`dim_A`, `dim_B`, `tol`, and per-batch `batch_<id>_attempted` /
`batch_<id>_accepted` counts). The generators store a `Vector{Matrix}` per
batch; `compare_detection.jl` stores a `Vector` of named tuples
`(state, robustness, min_dot, min_amp, dot_idx, amp_idx)` instead.

## Resumability and reproducibility

- **Resumable:** rerun a generator with the same `--output` to fill in only the
  missing batches; completed batches are detected and skipped.
- **Reproducible:** every candidate is seeded deterministically, and a batch
  keeps the lowest-index successes, so a given configuration produces the same
  dataset regardless of thread count. Changing `--batch` changes the seed layout
  and therefore the dataset.
