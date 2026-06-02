#!/usr/bin/env bash
## 4-chain parallel tMsPGOcc 500sp v3
## Memory-aware: starts with 2 chains to test, then expands
set -euo pipefail

PROJECT_ROOT="${BIRD_PROJECT_ROOT:-$HOME/bird_dynamic_occupancy_analysis}"
SCRIPT="$PROJECT_ROOT/code_v3/04_run_tMsPGOcc_500sp_v3.R"
LOG_DIR="$PROJECT_ROOT/logs_v3"
mkdir -p "$LOG_DIR"

N_CHAINS=4
N_OMP=32  # Reduced from 48 to save memory (32*4=128 cores, leaves 128 for system)
TS=$(date +%Y%m%d_%H%M%S)

echo "============================================================"
echo "[$(date)] Starting 4-chain tMsPGOcc 500sp v3"
echo "  Script: $SCRIPT"
echo "  OMP/chain: $N_OMP | Total OMP: $((N_OMP * N_CHAINS))"
echo "  RUN_LABEL: v3_full_500sp_ar1_temporal"
echo "  Memory check:"
free -h | head -3
echo "============================================================"

# Memory check: need at least 200GB free for 500sp
AVAIL_GB=$(free -g | awk "/Mem:/{print \$7}")
echo "[MEM] Available: ${AVAIL_GB}GB"
if [ "$AVAIL_GB" -lt 200 ]; then
  echo "[WARN] Only ${AVAIL_GB}GB available. 500sp model may need 300-500GB."
  echo "[WARN] Consider waiting for other processes to finish."
  echo "[WARN] Proceeding anyway - OOM possible."
fi

PIDS=()
for cid in 1 2 3 4; do
  cpu_start=$(( (cid - 1) * N_OMP ))
  cpu_end=$(( cpu_start + N_OMP - 1 ))
  cpu_range="${cpu_start}-${cpu_end}"
  log_file="$LOG_DIR/04_t500v3_chain${cid}_${TS}.log"
  echo "  [Chain $cid] CPU $cpu_range | Log: $log_file"
  BIRD_PROJECT_ROOT="$PROJECT_ROOT" \
  V3_CHAIN_ID=$cid \
  V3_N_OMP=$N_OMP \
  OMP_NUM_THREADS=$N_OMP \
  V3_MAX_SPECIES=500 \
  V3_RUN_LABEL=v3_full_500sp_ar1_temporal \
  taskset -c "$cpu_range" \
  Rscript "$SCRIPT" > "$log_file" 2>&1 &
  PIDS+=($!)
  # Small delay to stagger startup and reduce peak memory
  sleep 30
done

echo "  Waiting for ${#PIDS[@]} chain(s)... PIDs: ${PIDS[*]}"
echo "  Monitor: tail -f $LOG_DIR/04_t500v3_chain1_${TS}.log"

FAILED=0
for pid in "${PIDS[@]}"; do
  if ! wait "$pid"; then
    FAILED=$((FAILED + 1))
  fi
done
echo "[$(date)] Done. Failed: $FAILED/${#PIDS[@]}"

# After chains complete, combine them
if [ $FAILED -eq 0 ]; then
  echo "[$(date)] All chains succeeded. Running 04c_combine_chains..."
  cd "$PROJECT_ROOT"
  BIRD_PROJECT_ROOT="$PROJECT_ROOT" \
  V3_RUN_LABEL=v3_full_500sp_ar1_temporal \
  Rscript code_v3/04c_combine_chains.R > "$LOG_DIR/04c_t500v3_combine_${TS}.log" 2>&1
  echo "[$(date)] Chain combination done."
else
  echo "[ERROR] $FAILED chain(s) failed. Check logs before combining."
fi
