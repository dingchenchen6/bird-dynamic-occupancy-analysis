#!/usr/bin/env bash
## run_pipeline_05_to_13.sh -- 顺序执行 04b→05→06→06b→07→08→10→11→12→13
## 用法: nohup bash run_pipeline_05_to_13.sh > logs_v3/pipeline_05_13.log 2>&1 &

set -euo pipefail

PROJECT_ROOT="/home/dingchenchen/bird_dynamic_occupancy_analysis"
export BIRD_PROJECT_ROOT="$PROJECT_ROOT"
export V3_CODE_DIR="$PROJECT_ROOT/code_v3"

LOG_DIR="$PROJECT_ROOT/logs_v3"
MARKER_DIR="$LOG_DIR/markers"
mkdir -p "$LOG_DIR" "$MARKER_DIR"

RSCRIPT="$(command -v Rscript || echo /usr/local/bin/Rscript)"

# 检查04b是否已完成
wait_for_04b() {
  local max_wait=3600  # 最多等1小时
  local waited=0
  while true; do
    if [ -f "$MARKER_DIR/04b_done.marker" ]; then
      echo "[pipeline] 04b completed. Proceeding..."
      return 0
    fi
    if [ -f "$MARKER_DIR/04b_failed.marker" ]; then
      echo "[pipeline] 04b FAILED. Stopping pipeline."
      exit 1
    fi
    if [ $waited -ge $max_wait ]; then
      echo "[pipeline] Timeout waiting for 04b. Stopping."
      exit 1
    fi
    echo "[pipeline] Waiting for 04b... (${waited}s elapsed)"
    sleep 60
    waited=$((waited + 60))
  done
}

run_stage() {
  local stage_id="$1"
  local script_name="$2"
  local desc="$3"
  local log_file="${LOG_DIR}/${stage_id}_${script_name%.R}_$(date +%Y%m%d_%H%M%S).log"

  echo ""
  echo "============================================================"
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Stage ${stage_id} :: ${desc}"
  echo "  Script: ${script_name}"
  echo "  Log:    ${log_file}"
  echo "============================================================"

  local start_time=$(date +%s)
  if "$RSCRIPT" "${V3_CODE_DIR}/${script_name}" >"$log_file" 2>&1; then
    local end_time=$(date +%s)
    local elapsed=$((end_time - start_time))
    echo "  DONE in ${elapsed}s"
    echo "$(date '+%Y-%m-%d %H:%M:%S') ${stage_id} done" > "$MARKER_DIR/${stage_id}_done.marker"
  else
    local end_time=$(date +%s)
    local elapsed=$((end_time - start_time))
    echo "  FAILED after ${elapsed}s -- check ${log_file}" >&2
    echo "$(date '+%Y-%m-%d %H:%M:%S') ${stage_id} failed" > "$MARKER_DIR/${stage_id}_failed.marker"
    echo "  [WARNING] Continuing to next stage despite failure"
  fi
}

# 如果04b还没完成，等待
if [ ! -f "$MARKER_DIR/04b_done.marker" ] && [ ! -f "$MARKER_DIR/04b_failed.marker" ]; then
  wait_for_04b
fi

# 顺序执行后续阶段
run_stage "05" "05_postprocess_diversity.R" "Diversity post-processing (richness, PD, FDiv, FEve, Baselga)"
run_stage "06" "06_figures_publication.R" "Publication figures"
run_stage "06b" "06b_regenerate_driver_plots.R" "Driver plots regeneration"
run_stage "07" "07_render_manuscript.R" "Manuscript rendering"
run_stage "08" "08_post_full_run_pipeline.sh" "Post full-run pipeline"
run_stage "10" "10_finalize_with_dashline_pptx.R" "PPTX finalization"
run_stage "11" "11_render_journal_docx.R" "Journal docx rendering"
run_stage "12" "12_homogenization_spatiotemporal_maps.R" "Homogenization maps"
run_stage "13" "13_render_global_proposal_docx.R" "Global proposal docx"

echo ""
echo "============================================================"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Pipeline 05-13 COMPLETE"
echo "============================================================"
