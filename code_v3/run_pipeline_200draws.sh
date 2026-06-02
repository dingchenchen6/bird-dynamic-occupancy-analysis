#!/usr/bin/env bash
## run_pipeline_200draws.sh — 200 draws 版本 06-13 管线
## 与 400 draws 版本完全隔离，使用不同的 RUN_LABEL

set -uo pipefail

RUN_LABEL_200="v3_full_200sp_ar1_spatial_200draws"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" \u0026\u0026 pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
LOG_DIR="$PROJECT_ROOT/logs_v3"
mkdir -p "$LOG_DIR"

TS="$(date +%Y%m%d_%H%M%S)"
MASTER_LOG="$LOG_DIR/pipeline_200draws_${TS}.log"

RSCRIPT="$(which Rscript 2>/dev/null || echo /opt/R/4.5.3/lib64/R/bin/Rscript)"

# 设置环境变量，确保所有脚本使用 _200draws 后缀
export V3_RUN_LABEL="$RUN_LABEL_200"
export V3_CODE_DIR="$SCRIPT_DIR"
export BIRD_PROJECT_ROOT="$PROJECT_ROOT"

echo "========================================" | tee -a "$MASTER_LOG"
echo " 200 draws Pipeline: 06 → 13" | tee -a "$MASTER_LOG"
echo " RUN_LABEL: $RUN_LABEL_200" | tee -a "$MASTER_LOG"
echo " Started: $(date)" | tee -a "$MASTER_LOG"
echo "========================================" | tee -a "$MASTER_LOG"

FAILED=()

run_stage() {
  local num="$1"
  local script="$2"
  local desc="$3"
  local log="$LOG_DIR/${num}_200draws_${script%.R}_${TS}.log"

  echo "" | tee -a "$MASTER_LOG"
  echo "Stage $num: $desc" | tee -a "$MASTER_LOG"

  if "$RSCRIPT" "$SCRIPT_DIR/$script" > "$log" 2>\u00261; then
    echo "  ✅ Done" | tee -a "$MASTER_LOG"
  else
    echo "  ❌ FAILED" | tee -a "$MASTER_LOG"
    FAILED+=("$num")
  fi
}

# 06-13 顺序运行
run_stage "06" "06_figures_publication.R" "Publication figures"
run_stage "06b" "06b_regenerate_driver_plots.R" "Driver plots"
run_stage "07" "07_render_manuscript.R" "Manuscript"
run_stage "10" "10_finalize_with_dashline_pptx.R" "PPTX"
run_stage "11" "11_render_journal_docx.R" "Journal DOCX"
run_stage "12" "12_homogenization_spatiotemporal_maps.R" "Homogenization maps"
run_stage "13" "13_render_global_proposal_docx.R" "Global proposal"

echo "" | tee -a "$MASTER_LOG"
echo "========================================" | tee -a "$MASTER_LOG"
if [ ${#FAILED[@]} -eq 0 ]; then
  echo " 200 draws Pipeline ALL COMPLETE" | tee -a "$MASTER_LOG"
else
  echo " Failed stages: ${FAILED[*]}" | tee -a "$MASTER_LOG"
fi
echo " Finished: $(date)" | tee -a "$MASTER_LOG"
echo "========================================" | tee -a "$MASTER_LOG"
