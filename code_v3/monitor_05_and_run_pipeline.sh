#!/usr/bin/env bash
## monitor_05_and_run_pipeline.sh - Wait for PID to finish, then run 06-13

set -uo pipefail

PROJECT_ROOT="${BIRD_PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
LOG_DIR="$PROJECT_ROOT/logs_v3"
mkdir -p "$LOG_DIR"

PID_TO_WATCH=${1:-954637}
CHECK_INTERVAL=300
MASTER_LOG="$LOG_DIR/monitor_05_$(date +%Y%m%d_%H%M%S).log"

echo "======================================================" | tee "$MASTER_LOG"
echo " Monitoring PID $PID_TO_WATCH (05 postprocess)" | tee -a "$MASTER_LOG"
echo " Started: $(date)" | tee -a "$MASTER_LOG"
echo "======================================================" | tee -a "$MASTER_LOG"

while kill -0 "$PID_TO_WATCH" 2>/dev/null; do
  ELAPSED=$(ps -p "$PID_TO_WATCH" -o etimes= 2>/dev/null | tr -d " ")
  if [[ -n "$ELAPSED" ]]; then
    HOURS=$((ELAPSED / 3600))
    MINS=$(( (ELAPSED % 3600) / 60 ))
    echo "$(date +%H:%M:%S) PID $PID_TO_WATCH still running (${HOURS}h${MINS}m)" | tee -a "$MASTER_LOG"
  fi
  sleep $CHECK_INTERVAL
done

echo "" | tee -a "$MASTER_LOG"
echo "PID $PID_TO_WATCH finished at $(date)" | tee -a "$MASTER_LOG"

RUN_LABEL="v3_full_200sp_ar1_spatial"
PSI_THINNED="$PROJECT_ROOT/data/derived_v3/psi_samples_thinned_${RUN_LABEL}.rds"
DIVERSITY_SUMMARY="$PROJECT_ROOT/results_v3/table_diversity_summary_${RUN_LABEL}.csv"

if [[ -f "$PSI_THINNED" ]] && [[ -f "$DIVERSITY_SUMMARY" ]]; then
  echo "OK: 05 output files found" | tee -a "$MASTER_LOG"
else
  echo "WARN: Some 05 output files missing" | tee -a "$MASTER_LOG"
fi

echo "" | tee -a "$MASTER_LOG"
echo "Starting 06-13 pipeline" | tee -a "$MASTER_LOG"

export V3_RUN_LABEL="$RUN_LABEL"
export BIRD_PROJECT_ROOT="$PROJECT_ROOT"
export V3_CODE_DIR="$PROJECT_ROOT/code_v3"

bash "$PROJECT_ROOT/code_v3/98_smart_pipeline_05to13.sh" 2>&1 | tee -a "$MASTER_LOG"

echo "" | tee -a "$MASTER_LOG"
echo "======================================================" | tee -a "$MASTER_LOG"
echo " Monitor + Pipeline COMPLETE" | tee -a "$MASTER_LOG"
echo " Finished: $(date)" | tee -a "$MASTER_LOG"
echo "======================================================" | tee -a "$MASTER_LOG"
