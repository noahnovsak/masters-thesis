#!/usr/bin/env bash
# Assemble SUMMARY.md from the persisted run artifacts (logs, timings.tsv, ledgers).
# Idempotent: safe to re-run; reads only on-disk results. Run by the completion
# watcher when the driver prints "FINAL SCAN complete", but also fine to run early
# for a partial summary.
set -u
cd "$(dirname "$0")/../.."          # -> code/
RES=results/final_seed42
L=$RES/logs
OUT=$RES/SUMMARY.md

sect() { grep -aA"$2" -- "$3" "$1" 2>/dev/null | sed 's/\x1b\[[0-9;]*m//g'; }   # <file> <lines> <pattern>; strip ANSI
line() { grep -a "$2" "$1" 2>/dev/null | sed 's/\x1b\[[0-9;]*m//g' | tail -1; }
ledstats() { local f=$1; [ -f "$f" ] || { echo "n/a"; return; }
  awk -F, 'NR>1{t++; if($3=="true")d++} END{printf "%d tested, %d detected", t+0, d+0}' "$f"; }

{
echo "# PPT² final scan (seed 42, 4×4) — results summary"
echo
echo "_Generated $(date '+%F %T') from on-disk artifacts in \`$RES/\`._"
echo "_Methodology and prose: see [WRITEUP.md](WRITEUP.md)._"
echo
grep -qa "FINAL SCAN complete" "$RES/driver.out" 2>/dev/null \
  && echo "**Status: COMPLETE.**" || echo "**Status: PARTIAL (run still in progress).**"
echo
echo "## Configuration"
echo "- n = m = 4, DPS level 2, tol = 1e-8, base seed 42, -t 40 (Mosek/BLAS defaults)."
echo "- Random-PPT pools capped at 5000 states each; PnCP library, witness states, witness-PPT² and pair tests at full 10000."
echo
echo "## Timings"
echo '```'
column -t -s$'\t' "$RES/timings.tsv" 2>/dev/null
echo '```'
echo
echo "## 1. Libraries"
echo "**PnCP witness maps** (\`gen_pncp\`): $(line "$L/01_gen_pncp.log" 'Done. Accepted')"
echo
echo "**Witness-constructed PPT states** (\`gen_witness_ppt\`):"
echo '```'; sect "$L/02_gen_witness_ppt.log" 2 'Witness-restricted PPT minimisation'; echo '```'
echo
echo "## 2. Detection power — random PPT pools (5000 each)"
echo "### Asymmetric"
echo '```'; sect "$L/03_cmp_asym.log" 8 'Detection efficacy over'; echo '```'
echo "### Symmetric"
echo '```'; sect "$L/04_cmp_sym.log" 8 'Detection efficacy over'; echo '```'
echo
echo "## 3. Detection power — witness-constructed pool (witness-vs-DPS comparison)"
echo '```'; sect "$L/05_detection_power_witness.log" 11 'Detection power over'; echo '```'
echo
echo "## 4. PPT² conjecture via the witness SDP (all 10000 witnesses)"
echo '```'; sect "$L/06_gen_witness_ppt2.log" 4 'See-saw PPT² witness search'; echo '```'
echo
echo "## 5. PPT² on 10000 compositions per state family (three criteria + DPS)"
echo "| pool | result | ledger |"
echo "|---|---|---|"
echo "| witness    | $(line "$L/07_test_ppt2_witness.log" 'Cumulative:') | $(ledstats "$RES/ppt2_results/witness/tested_witness_ppt_4x4.csv") |"
echo "| asymmetric | $(line "$L/08_test_ppt2_asym.log"    'Cumulative:') | $(ledstats "$RES/ppt2_results/asym/tested_cmp_asym_4x4.csv") |"
echo "| symmetric  | $(line "$L/09_test_ppt2_sym.log"     'Cumulative:') | $(ledstats "$RES/ppt2_results/sym/tested_cmp_sym_4x4.csv") |"
echo
echo "Any detected composition (a PPT² counterexample) is saved as \`ppt2_results/<pool>/result_<i>_<j>.jld2\`."
} > "$OUT"
echo "wrote $OUT"
