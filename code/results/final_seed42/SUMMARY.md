# PPT² final scan (seed 42, 4×4) — results summary

_Generated 2026-06-08 09:37:25 from on-disk artifacts in `results/final_seed42/`._
_Methodology and prose: see [WRITEUP.md](WRITEUP.md)._

**Status: COMPLETE.**

## Configuration
- n = m = 4, DPS level 2, tol = 1e-8, base seed 42, -t 40 (Mosek/BLAS defaults).
- Random-PPT pools capped at 5000 states each; PnCP library, witness states, witness-PPT² and pair tests at full 10000.

## Timings
```
step                     seconds  hh:mm:ss  start                end                  note1       note2
gen_pncp                 2352     00:39:12  2026-06-06 11:01:56  2026-06-06 11:41:08  attempts=1  rc=0
gen_witness_ppt          214      00:03:34  2026-06-06 11:41:08  2026-06-06 11:44:42  attempts=1  rc=0
compare_detection_asym   51       00:00:51  2026-06-07 09:02:16  2026-06-07 09:03:07  attempts=1  rc=0
compare_detection_sym    16264    04:31:04  2026-06-07 09:03:07  2026-06-07 13:34:11  attempts=1  rc=0
detection_power_witness  21896    06:04:56  2026-06-07 13:34:11  2026-06-07 19:39:07  attempts=1  rc=0
gen_witness_ppt2         5610     01:33:30  2026-06-07 19:39:07  2026-06-07 21:12:37  attempts=1  rc=0
test_ppt2_witness        12735    03:32:15  2026-06-07 21:12:37  2026-06-08 00:44:52  attempts=1  rc=0
test_ppt2_asym           16599    04:36:39  2026-06-08 00:44:52  2026-06-08 05:21:31  attempts=1  rc=0
test_ppt2_sym            15334    04:15:34  2026-06-08 05:21:31  2026-06-08 09:37:05  attempts=1  rc=0
```

## 1. Libraries
**PnCP witness maps** (`gen_pncp`): Done. Accepted 10000/10000 (100.0%).

**Witness-constructed PPT states** (`gen_witness_ppt`):
```
Witness-restricted PPT minimisation over 10000 witnesses (tol 1.0e-8):
  certified a PPT entangled state: 10000/10000 (100.0%)
  optimum tr(W·ρ):  min -0.03484, median -5.292e-5, max -6.057e-8
```

## 2. Detection power — random PPT pools (5000 each)
### Asymmetric
```
Detection efficacy over 5000 states (each detected by ≥1 criterion):
  DPS (level 2)               5000   (100.0%)
  PnCP trace witness          0      (0.0%)
  PnCP ampliation             0      (0.0%)
  PnCP (trace OR ampliation)  0      (0.0%)
  ------------------------------------
  DPS only (PnCP missed)      5000   (100.0%)
  PnCP only (DPS missed)      0      (0.0%)
  DPS and PnCP                0      (0.0%)
```
### Symmetric
```
Detection efficacy over 5000 states (each detected by ≥1 criterion):
  DPS (level 2)               5000   (100.0%)
  PnCP trace witness          0      (0.0%)
  PnCP ampliation             0      (0.0%)
  PnCP (trace OR ampliation)  0      (0.0%)
  ------------------------------------
  DPS only (PnCP missed)      5000   (100.0%)
  PnCP only (DPS missed)      0      (0.0%)
  DPS and PnCP                0      (0.0%)
```

## 3. Detection power — witness-constructed pool (witness-vs-DPS comparison)
```
Detection power over 10000 states from witness_ppt_4x4.jld2:
  DPS (level 2)               10000  (100.0%)
  PnCP trace witness          10000  (100.0%)
  PnCP ampliation (system=1)  10000  (100.0%)
  PnCP (trace OR ampliation)  10000  (100.0%)
  ANY criterion               10000  (100.0%)
  ------------------------------------
  DPS only (PnCP missed)      0      (0.0%)
  PnCP only (DPS missed)      0      (0.0%)
  DPS and PnCP                10000  (100.0%)
```

## 4. PPT² conjecture via the witness SDP (all 10000 witnesses)
```
See-saw PPT² witness search over 10000 witnesses (tol 1.0e-8):
  detected a composition (PPT² counterexample candidate): 0/10000 (0.0%)
  see-saw optimum tr(W·composite):  min 1.557e-16, median 2.985e-10, max 2.237e-9
  (PPT² predicts every optimum ≥ 0; a value below -tol is a counterexample candidate.)
```

## 5. PPT² on 10000 compositions per state family (three criteria + DPS)
| pool | result | ledger |
|---|---|---|
| witness    | Cumulative: 10000 tested, 0 detected (of 10000 possible pairs). | 10000 tested, 0 detected |
| asymmetric | Cumulative: 10000 tested, 0 detected (of 10000 possible pairs). | 10000 tested, 0 detected |
| symmetric  | Cumulative: 10000 tested, 0 detected (of 10000 possible pairs). | 10000 tested, 0 detected |

Any detected composition (a PPT² counterexample) is saved as `ppt2_results/<pool>/result_<i>_<j>.jld2`.
