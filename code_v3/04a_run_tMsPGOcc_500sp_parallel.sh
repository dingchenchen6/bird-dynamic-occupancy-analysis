#!/usr/bin/env bash
## 4-chain parallel tMsPGOcc 500sp
set -euo pipefail
PROJECT_ROOT="${BIRD_PROJECT_ROOT:-$HOME/bird_dynamic_occupancy_analysis}"
SCRIPT="$PROJECT_ROOT/code_v3/04_run_tMsPGOcc_500sp_v2.R"
LOG_DIR="$PROJECT_ROOT/logs_v3"
mkdir -p "$LOG_DIR"

N_CHAINS=4
N_OMP=48
TS=$(date +%Y%m%d_%H%M%S)
echo "[$(date)] Starting 4-chain tMsPGOcc 500sp, OMP/chain=$N_OMP"

PIDS=()
for cid in 1 2 3 4; do
  cpu_start=$(( (cid - 1) * N_OMP ))
  cpu_end=$(( cpu_start + N_OMP - 1 ))
  cpu_range="${cpu_start}-${cpu_end}"
  log_file="$LOG_DIR/04_t500_chain${cid}_${TS}.log"
  echo "  [Chain $cid] CPU $cpu_range | Log: $log_file"
  V3_CHAIN_ID=$cid V3_N_OMP=$N_OMP V3_MAX_SPECIES=500 V3_RUN_LABEL=v3_full_500sp_ar1_temporal \
    OMP_NUM_THREADS=$N_OMP \
    taskset -c "$cpu_range" \
    Rscript "$SCRIPT" > "$log_file" 2>&1 &
  PIDS+=($!)
done

echo "  Waiting for ${#PIDS[@]} chain(s)... PIDs: ${PIDS[*]}"
FAILED=0
for pid in "${PIDS[@]}"; do
  if ! wait "$pid"; then
    FAILED=$((FAILED + 1))
  fi
done
echo "[$(date)] Done. Failed: $FAILED/${#PIDS[@]}"
