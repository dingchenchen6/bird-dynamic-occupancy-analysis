#!/usr/bin/env bash
## auto_monitor_and_pipeline.sh -- 监控 04b 并在完成后自动顺序执行 05-13

set -euo pipefail

PROJECT_ROOT="/home/dingchenchen/bird_dynamic_occupancy_analysis"
export BIRD_PROJECT_ROOT="$PROJECT_ROOT"
export V3_CODE_DIR="$PROJECT_ROOT/code_v3"

LOG_DIR="$PROJECT_ROOT/logs_v3"
MARKER_DIR="$LOG_DIR/markers"
mkdir -p "$MARKER_DIR"

RSCRIPT="$(command -v Rscript || echo /usr/local/bin/Rscript)"
MONITOR_LOG="$LOG_DIR/auto_monitor_$(date +%Y%m%d_%H%M%S).log"

log_msg() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$MONITOR_LOG"
}

# 等待 04b 完成（检查合并后的 fit 文件或标记）
wait_for_04b() {
  local fit_file="$PROJECT_ROOT/data/derived_v3/tMsPGOcc_fit_v3_full_200sp_ar1_spatial.rds"
  local max_wait_hours=12
  local waited=0

  log_msg "Waiting for 04b to complete (max ${max_wait_hours}h)..."

  while true; do
    # 检查合并后的 fit 文件是否存在且大小合理
    if [[ -f "$fit_file" ]]; then
      local fsize=$(stat -c%s "$fit_file" 2>/dev/null || echo 0)
      if [[ "$fsize" -gt 1073741824 ]]; then  # > 1GB
        log_msg "04b complete: $fit_file ($(numfmt --to=iec $fsize))"
        echo "$(date '+%Y-%m-%d %H:%M:%S') 04b_done" > "$MARKER_DIR/04b_done.marker"
        return 0
      fi
    fi

    # 检查是否有 04b 的 done marker
    if [[ -f "$MARKER_DIR/04b_done.marker" ]]; then
      log_msg "04b marker found."
      return 0
    fi

    # 检查是否超时
    if [[ $waited -ge $((max_wait_hours * 3600)) ]]; then
      log_msg "TIMEOUT waiting for 04b after ${max_wait_hours}h"
      return 1
    fi

    log_msg "Still waiting for 04b... (${waited}s elapsed)"
    sleep 300  # 5分钟检查一次
    waited=$((waited + 300))
  done
}

run_stage() {
  local stage_id="$1"
  local script_name="$2"
  local desc="$3"
  local log_file="${LOG_DIR}/${stage_id}_${script_name%.R}_$(date +%Y%m%d_%H%M%S).log"

  log_msg "============================================================"
  log_msg "Stage ${stage_id} :: ${desc}"
  log_msg "  Script: ${script_name}"
  log_msg "  Log:    ${log_file}"
  log_msg "============================================================"

  local start_time=$(date +%s)
  if "$RSCRIPT" "${V3_CODE_DIR}/${script_name}" >"$log_file" 2>&1; then
    local end_time=$(date +%s)
    local elapsed=$((end_time - start_time))
    log_msg "Stage ${stage_id} DONE in ${elapsed}s"
    echo "$(date '+%Y-%m-%d %H:%M:%S') ${stage_id} done" > "$MARKER_DIR/${stage_id}_done.marker"
  else
    local end_time=$(date +%s)
    local elapsed=$((end_time - start_time))
    log_msg "Stage ${stage_id} FAILED after ${elapsed}s -- check ${log_file}"
    echo "$(date '+%Y-%m-%d %H:%M:%S') ${stage_id} failed" > "$MARKER_DIR/${stage_id}_failed.marker"
    log_msg "Continuing to next stage..."
  fi
}

# 主流程
log_msg "Auto-monitor started."

if wait_for_04b; then
  # 顺序执行后续阶段
  run_stage "05" "05_postprocess_diversity.R" "Diversity post-processing"
  run_stage "06" "06_figures_publication.R" "Publication figures"
  run_stage "06b" "06b_regenerate_driver_plots.R" "Driver plots"
  run_stage "07" "07_render_manuscript.R" "Manuscript rendering"
  run_stage "08" "08_post_full_run_pipeline.sh" "Post full-run pipeline"
  run_stage "10" "10_finalize_with_dashline_pptx.R" "PPTX finalization"
  run_stage "11" "11_render_journal_docx.R" "Journal docx"
  run_stage "12" "12_homogenization_spatiotemporal_maps.R" "Homogenization maps"
  run_stage "13" "13_render_global_proposal_docx.R" "Global proposal docx"

  log_msg "============================================================"
  log_msg "ALL STAGES COMPLETE"
  log_msg "============================================================"
else
  log_msg "Pipeline aborted: 04b did not complete."
  exit 1
fi
