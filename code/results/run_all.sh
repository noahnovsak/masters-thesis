#!/usr/bin/env bash
# Final-scan driver: PPT2 thesis pipeline, seed 42, 4x4. Steps run sequentially so
# each Julia process gets the full thread budget (oversubscribing OOMs on the level-2
# DPS solves). Every step logs to logs/, appends its wall time to timings.tsv, is
# retried 3x, and is marked done with a .done_<step> file so a rerun resumes.
set -u

cd "$(dirname "$0")/.."
RES=results
LOGS=$RES/logs
THREADS=40
# Resolve the Julia binary directly to avoid juliaup lock
JBIN=$(ls -d "$HOME"/.julia/juliaup/julia-*/bin/julia 2>/dev/null | head -1)
JBIN=${JBIN:-julia}
J="$JBIN --project=. -t $THREADS"

FORMS=$RES/pncp_4x4.jld2
ASYM=$RES/cmp_asym_4x4.jld2
SYM=$RES/cmp_sym_4x4.jld2
WPPT=$RES/witness_ppt_4x4.jld2
WPPT2=$RES/witness_ppt2_4x4.jld2
DP_WIT=$RES/detection_power_witness_4x4.jld2

TIMINGS=$RES/timings.tsv
[ -f "$TIMINGS" ] || printf "step\tseconds\thh:mm:ss\tstart\tend\n" > "$TIMINGS"

run_step () {                        # run_step <name> <logfile> <cmd...>
    local name=$1 log=$2; shift 2
    local attempts=3 a rc start end t0 t1 dur total=0
    if [ -f "$RES/.done_$name" ]; then
        echo ""
        echo ">>> [$name] already complete (marker $RES/.done_$name) — skipping."
        return 0
    fi
    start=$(date '+%F %T')
    echo ""
    echo "================================================================"
    echo ">>> [$name] starting $(date '+%F %T')"
    echo ">>> cmd: $*"
    echo "================================================================"
    for ((a=1; a<=attempts; a++)); do
        [ $a -gt 1 ] && echo ">>> [$name] retry $a/$attempts (resumable steps continue from checkpoint) $(date '+%F %T')"
        t0=$(date +%s)
        "$@" > "$log" 2>&1
        rc=$?
        t1=$(date +%s); dur=$((t1 - t0)); total=$((total + dur))
        [ $rc -eq 0 ] && break
        echo "!!! [$name] attempt $a/$attempts failed (rc=$rc) after ${dur}s"
        tail -15 "$log"
        sleep 10
    done
    end=$(date '+%F %T')
    printf "%s\t%d\t%02d:%02d:%02d\t%s\t%s\tattempts=%d\trc=%d\n" \
        "$name" "$total" $((total/3600)) $(((total%3600)/60)) $((total%60)) "$start" "$end" "$a" "$rc" >> "$TIMINGS"
    if [ $rc -ne 0 ]; then
        echo "!!! [$name] FAILED after $attempts attempts (rc=$rc, ${total}s total) — see $log"
        exit $rc
    fi
    touch "$RES/.done_$name"
    echo ">>> [$name] done in ${total}s (attempt $a). tail:"
    tail -6 "$log"
}

echo "##### FINAL SCAN (seed 42) started $(date '+%F %T') on $THREADS threads #####"

# 1. PnCP witness-map library (seed 42) — blocks everything downstream.
run_step gen_pncp "$LOGS/01_gen_pncp.log" \
    $J scripts/gen_pncp.jl --total 10000 --batch 200 -n 4 -m 4 --seed 42 -o "$FORMS"

# 2. Witness-constructed bound-entangled PPT states (one SDP per form; fast).
run_step gen_witness_ppt "$LOGS/02_gen_witness_ppt.log" \
    $J scripts/gen_witness_ppt.jl -n 4 -m 4 --tol 1e-8 -f "$FORMS" -o "$WPPT"

# 3. Random PPT pool, ASYMMETRIC: sample + record DPS/trace/ampliation detection power.
run_step compare_detection_asym "$LOGS/03_cmp_asym.log" \
    $J scripts/compare_detection.jl --total 5000 --batch 200 -n 4 -m 4 \
       --level 2 --tol 1e-8 --seed 42 -f "$FORMS" -o "$ASYM"

# 4. Random PPT pool, SYMMETRIC (partial-transpose invariant).
run_step compare_detection_sym "$LOGS/04_cmp_sym.log" \
    $J scripts/compare_detection.jl --total 5000 --batch 200 -n 4 -m 4 \
       --level 2 --tol 1e-8 --seed 42 --ppt-invariant -f "$FORMS" -o "$SYM"

# 5. Detection power of the three tests on the witness-constructed library.
run_step detection_power_witness "$LOGS/05_detection_power_witness.log" \
    $J scripts/detection_power.jl -n 4 -m 4 --level 2 --tol 1e-8 -s "$WPPT" -f "$FORMS" -o "$DP_WIT"

# 6. PPT2 conjecture via the witness SDP (see-saw over PPT-map compositions), all witnesses.
run_step gen_witness_ppt2 "$LOGS/06_gen_witness_ppt2.log" \
    $J scripts/gen_witness_ppt2.jl -n 4 -m 4 --tol 1e-8 --restarts 16 --max_iter 40 \
       --seed 42 -f "$FORMS" -o "$WPPT2"

# 7-9. PPT2 conjecture on 100x100 = 10000 ordered pairs, all three criteria (+dps),
#      for each of the three state pools. Separate output dirs so result files / ledgers
#      never collide.
run_step test_ppt2_witness "$LOGS/07_test_ppt2_witness.log" \
    $J scripts/test_ppt2.jl -n 4 -m 4 --tol 1e-8 -s "$WPPT" -f "$FORMS" \
       --max-states 100 --with-dps -o "$RES/ppt2_results/witness"

run_step test_ppt2_asym "$LOGS/08_test_ppt2_asym.log" \
    $J scripts/test_ppt2.jl -n 4 -m 4 --tol 1e-8 -s "$ASYM" -f "$FORMS" \
       --max-states 100 --with-dps -o "$RES/ppt2_results/asym"

run_step test_ppt2_sym "$LOGS/09_test_ppt2_sym.log" \
    $J scripts/test_ppt2.jl -n 4 -m 4 --tol 1e-8 -s "$SYM" -f "$FORMS" \
       --max-states 100 --with-dps -o "$RES/ppt2_results/sym"

# 10-11. Cross-detection: test every witness-derived state against the WHOLE witness
#        library (10^8 trace pairs, plus the stronger ampliation test). Both read the
#        pncp/witness libraries from DATADIR and write cross_*.jld2 there.
run_step cross_trace "$LOGS/10_cross_trace.log" \
    env DATADIR="$RES" $J scripts/cross_trace.jl

run_step cross_ampl "$LOGS/11_cross_ampl.log" \
    env DATADIR="$RES" $J scripts/cross_ampl.jl

echo ""
echo "##### FINAL SCAN complete $(date '+%F %T') #####"
cat "$TIMINGS"
