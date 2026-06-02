#!/usr/bin/env bash
## 04c_parallel_and_full_pipeline.sh - 04c model + 05 postprocess + 06-13

set -euo pipefail

PROJECT_ROOT="${BIRD_PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
SCRIPT_DIR="$PROJECT_ROOT/code_v3"
LOG_DIR="$PROJECT_ROOT/logs_v3"
mkdir -p "$LOG_DIR"

RUN_LABEL_04C="v3_full_200sp_ar1_spatial_04c"
N_CHAINS=4

TOTAL_CORES=$(nproc 2>/dev/null || echo 256)
RESERVED=4
N_OMP_PER_CHAIN=$(( (TOTAL_CORES - RESERVED) / N_CHAINS ))
N_OMP_PER_CHAIN=$(( N_OMP_PER_CHAIN < 8 ? 8 : N_OMP_PER_CHAIN ))
N_OMP_PER_CHAIN=$(( N_OMP_PER_CHAIN > 63 ? 63 : N_OMP_PER_CHAIN ))

MASTER_LOG="$LOG_DIR/pipeline_04c_$(date +%Y%m%d_%H%M%S).log"

echo "======================================================" | tee "$MASTER_LOG"
echo " 04c + Full Pipeline" | tee -a "$MASTER_LOG"
echo " Run label: $RUN_LABEL_04C" | tee -a "$MASTER_LOG"
echo " Total cores: $TOTAL_CORES, OMP/chain: $N_OMP_PER_CHAIN" | tee -a "$MASTER_LOG"
echo " Started: $(date)" | tee -a "$MASTER_LOG"
echo "======================================================" | tee -a "$MASTER_LOG"

# Phase 1: 04c 4-chain parallel
echo "" | tee -a "$MASTER_LOG"
echo "-- Phase 1: 04c 4-chain parallel --" | tee -a "$MASTER_LOG"

SKIP_CHAINS=()
for CHAIN_ID in $(seq 1 $N_CHAINS); do
  CHAIN_FILE="$PROJECT_ROOT/data/derived_v3/tMsPGOcc_fit_${RUN_LABEL_04C}_chain${CHAIN_ID}.rds"
  if [[ -f "$CHAIN_FILE" ]]; then
    FSIZE=$(stat -c%s "$CHAIN_FILE" 2>/dev/null || echo 0)
    if [[ "$FSIZE" -gt 10485760 ]]; then
      echo "  Chain $CHAIN_ID already exists - skipping" | tee -a "$MASTER_LOG"
      SKIP_CHAINS+=("$CHAIN_ID")
    else
      rm -f "$CHAIN_FILE"
    fi
  fi
done

PIDS=()
CHAIN_IDS=()

for CHAIN_ID in $(seq 1 $N_CHAINS); do
  if [[ " ${SKIP_CHAINS[*]} " == *" $CHAIN_ID "* ]]; then continue; fi

  CPU_START=$(( (CHAIN_ID - 1) * N_OMP_PER_CHAIN ))
  CPU_END=$(( CPU_START + N_OMP_PER_CHAIN - 1 ))
  CPU_END=$(( CPU_END >= TOTAL_CORES ? TOTAL_CORES - 1 : CPU_END ))
  AFFINITY="${CPU_START}-${CPU_END}"

  CHAIN_LOG="$LOG_DIR/04c_chain${CHAIN_ID}_$(date +%Y%m%d_%H%M%S).log"
  echo "  [Chain $CHAIN_ID] CPU: $AFFINITY, OMP=$N_OMP_PER_CHAIN" | tee -a "$MASTER_LOG"

  V3_CHAIN_ID="$CHAIN_ID" \
  V3_CODE_DIR="$SCRIPT_DIR" \
  BIRD_PROJECT_ROOT="$PROJECT_ROOT" \
  V3_N_OMP_THREADS="$N_OMP_PER_CHAIN" \
  OMP_NUM_THREADS="$N_OMP_PER_CHAIN" \
  OMP_PROC_BIND=close \
  OMP_PLACES=threads \
  taskset -c "$AFFINITY" \
    Rscript "${SCRIPT_DIR}/04c_run_stMsPGOcc_4chain_merge.R" >"$CHAIN_LOG" 2>&1 &

  PIDS+=($!)
  CHAIN_IDS+=("$CHAIN_ID")
done

if [[ ${#PIDS[@]} -gt 0 ]]; then
  echo "  Waiting for ${#PIDS[@]} chain(s)..." | tee -a "$MASTER_LOG"
  FAILED=0
  for i in "${!PIDS[@]}"; do
    PID="${PIDS[$i]}"
    CHAIN="${CHAIN_IDS[$i]}"
    if wait "$PID"; then
      echo "  [Chain $CHAIN] Done" | tee -a "$MASTER_LOG"
    else
      echo "  [Chain $CHAIN] FAILED" | tee -a "$MASTER_LOG"
      FAILED=$((FAILED + 1))
    fi
  done

  if [[ $FAILED -gt 0 ]]; then
    echo "  $FAILED chain(s) failed! Aborting." | tee -a "$MASTER_LOG"
    exit 1
  fi
else
  echo "  All chains already completed." | tee -a "$MASTER_LOG"
fi

echo "  Phase 1 done: $(date)" | tee -a "$MASTER_LOG"

# Phase 2: 04c combine 4 chains
echo "" | tee -a "$MASTER_LOG"
echo "-- Phase 2: 04c combine 4 chains --" | tee -a "$MASTER_LOG"

COMBINE_LOG="$LOG_DIR/04c_combine_$(date +%Y%m%d_%H%M%S).log"
if Rscript "${SCRIPT_DIR}/04c_combine_chains.R" >"$COMBINE_LOG" 2>&1; then
  echo "  OK: Combine done" | tee -a "$MASTER_LOG"
else
  echo "  FAIL: Combine FAILED" | tee -a "$MASTER_LOG"
  exit 1
fi

echo "  Phase 2 done: $(date)" | tee -a "$MASTER_LOG"

# Phase 3: 05 postprocess (using 04c output)
echo "" | tee -a "$MASTER_LOG"
echo "-- Phase 3: 05 postprocess diversity (04c) --" | tee -a "$MASTER_LOG"

export V3_RUN_LABEL="$RUN_LABEL_04C"
export BIRD_PROJECT_ROOT="$PROJECT_ROOT"
export V3_CODE_DIR="$SCRIPT_DIR"

POST_LOG="$LOG_DIR/05_postprocess_${RUN_LABEL_04C}_$(date +%Y%m%d_%H%M%S).log"
if Rscript "${SCRIPT_DIR}/05_postprocess_diversity.R" >"$POST_LOG" 2>&1; then
  echo "  OK: 05 done" | tee -a "$MASTER_LOG"
else
  echo "  FAIL: 05 FAILED" | tee -a "$MASTER_LOG"
fi

echo "  Phase 3 done: $(date)" | tee -a "$MASTER_LOG"

# Phase 4: 06-13 publication pipeline
echo "" | tee -a "$MASTER_LOG"
echo "-- Phase 4: 06-13 publication pipeline (04c) --" | tee -a "$MASTER_LOG"

bash "${SCRIPT_DIR}/98_smart_pipeline_05to13.sh" 2>&1 | tee -a "$MASTER_LOG"

echo "" | tee -a "$MASTER_LOG"
echo "======================================================" | tee -a "$MASTER_LOG"
echo " 04c Full Pipeline COMPLETE" | tee -a "$MASTER_LOG"
echo " Finished: $(date)" | tee -a "$MASTER_LOG"
echo "======================================================" | tee -a "$MASTER_LOG"
