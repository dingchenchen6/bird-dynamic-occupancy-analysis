#!/usr/bin/env bash
## 91_run_v3_stable.sh  —  v3 稳定运行脚本（断网不中断）
##
## 改进:
##   1. 先 100km 再 10km（100km 可复用 v2 数据，10km 需重新提取环境变量）
##   2. 使用 nohup + tmux 双重保护，SSH 断开不影响
##   3. 每个 grid size 独立完成全流水线（01-08）
##   4. 模型保存失败自动重试 + checkpoint
##   5. 关键阶段前后写进度标记文件
##
## 用法:
##   # 方式 1: 在 tmux 内运行（推荐）
##   tmux new -s bird_v3
##   bash code_v3/91_run_v3_stable.sh
##
##   # 方式 2: nohup 后台运行
##   nohup bash code_v3/91_run_v3_stable.sh > logs_v3/pipeline_stable.log 2>&1 &

set -euo pipefail

# ── 配置 ───────────────────────────────────────────────────────────────
PROJECT_ROOT="/home/dingchenchen/bird_dynamic_occupancy_analysis"
export BIRD_PROJECT_ROOT="$PROJECT_ROOT"
export V3_CODE_DIR="$PROJECT_ROOT/code_v3"

LOG_DIR="$PROJECT_ROOT/logs_v3"
MARKER_DIR="$PROJECT_ROOT/logs_v3/markers"
mkdir -p "$LOG_DIR" "$MARKER_DIR" "$PROJECT_ROOT/data/derived_v3" \
  "$PROJECT_ROOT/results_v3" "$PROJECT_ROOT/figures_v3"

RSCRIPT="$(command -v Rscript || echo /usr/local/bin/Rscript)"
PIPELINE_LOG="$LOG_DIR/pipeline_stable_$(date +%Y%m%d_%H%M%S).log"

# 同时输出到屏幕和日志文件
exec > >(tee -a "$PIPELINE_LOG") 2>&1

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

write_marker() {
  local name="$1"
  echo "$(timestamp) $name" > "$MARKER_DIR/${name}.marker"
}

# ── 运行单个阶段 ─────────────────────────────────────────────────────
run_stage() {
  local stage_id="$1"
  local script_name="$2"
  local desc="$3"
  local extra_env="${4:-}"
  local log_file="${LOG_DIR}/${stage_id}_${script_name%.R}_$(date +%Y%m%d_%H%M%S).log"

  echo ""
  echo "============================================================"
  echo "[$(timestamp)] Stage ${stage_id} :: ${desc}"
  echo "  Script: ${script_name}"
  echo "  Log:    ${log_file}"
  echo "============================================================"

  write_marker "${stage_id}_start"

  local start_time=$(date +%s)
  if env $extra_env "$RSCRIPT" "${V3_CODE_DIR}/${script_name}" >"$log_file" 2>&1; then
    local end_time=$(date +%s)
    local elapsed=$(( end_time - start_time ))
    echo "  DONE in ${elapsed}s"
    write_marker "${stage_id}_done"
  else
    local end_time=$(date +%s)
    local elapsed=$(( end_time - start_time ))
    echo "  FAILED after ${elapsed}s — check ${log_file}" >&2
    write_marker "${stage_id}_failed"
    # 不退出，继续下一阶段（允许部分失败后手动恢复）
    echo "  [WARNING] Continuing to next stage despite failure"
  fi
}

# ── 运行单个 grid size 的完整流水线 ─────────────────────────────────
run_pipeline_for_grid() {
  local grid_km="$1"
  local run_label="$2"
  local is_pilot="$3"

  export V3_RUN_LABEL="$run_label"
  export V3_PILOT="$is_pilot"
  export V3_GRID_SIZE_KM="$grid_km"

  # 更新 RUN_LABEL 以包含网格尺寸
  if [[ "$grid_km" != "100" ]]; then
    export V3_RUN_LABEL="${run_label}_${grid_km}km"
  fi

  echo ""
  echo "############################################################"
  echo "[$(timestamp)] Starting pipeline for ${grid_km}km grid"
  echo "  Run label: $V3_RUN_LABEL"
  echo "  Pilot:     $is_pilot"
  echo "############################################################"

  # ── Stage 01-03: 前处理 ──────────────────────────────────────────
  run_stage "01_${grid_km}km" "01_merge_birdwatch_ebird.R" \
    "Merge and dedup observation records"

  run_stage "02_${grid_km}km" "02_build_survey_history.R" \
    "Build breeding-season survey history (${grid_km}km)"

  run_stage "03_${grid_km}km" "03_prepare_environment.R" \
    "Prepare environment (${grid_km}km)"

  run_stage "03b_${grid_km}km" "03b_extend_traits.R" \
    "Extend functional traits"

  # Stage 03c-03e: 环境变化数据（如果已存在则跳过）
  if [[ ! -f "$PROJECT_ROOT/data/derived_v3/climate_change_v3.rds" ]]; then
    run_stage "03c_${grid_km}km" "03c_prepare_climate_change.R" \
      "Prepare climate change covariates"
  else
    echo "[$(timestamp)] Skipping 03c — climate_change_v3.rds already exists"
  fi

  if [[ ! -f "$PROJECT_ROOT/data/derived_v3/landuse_change_v3.rds" ]]; then
    run_stage "03d_${grid_km}km" "03d_prepare_landuse_change.R" \
      "Prepare land use change covariates"
  else
    echo "[$(timestamp)] Skipping 03d — landuse_change_v3.rds already exists"
  fi

  if [[ ! -f "$PROJECT_ROOT/data/derived_v3/hfi_change_v3.rds" ]]; then
    run_stage "03e_${grid_km}km" "03e_prepare_hfi_change.R" \
      "Prepare HFI change covariates"
  else
    echo "[$(timestamp)] Skipping 03e — hfi_change_v3.rds already exists"
  fi

  # ── Stage 04: stMsPGOcc 模型（4链并行 + 降低OMP线程防segfault）───
  echo ""
  echo "============================================================"
  echo "[$(timestamp)] Stage 04_${grid_km}km :: stMsPGOcc 4-chain PARALLEL fit"
  echo "============================================================"
  write_marker "04_${grid_km}km_start"

  # 稳定性设置：降低每链OMP线程，绑定物理核心
  export V3_N_OMP_PER_CHAIN=12
  export V3_N_CHAINS=4
  export OMP_PROC_BIND=close
  export OMP_PLACES=cores

  local start_time=$(date +%s)
  if bash "${V3_CODE_DIR}/04a_run_stMsPGOcc_parallel.sh" \
      >"${LOG_DIR}/04_${grid_km}km_parallel_$(date +%Y%m%d_%H%M%S).log" 2>&1; then
    local end_time=$(date +%s)
    local elapsed=$(( end_time - start_time ))
    echo "  DONE in ${elapsed}s"
    write_marker "04_${grid_km}km_done"
  else
    local end_time=$(date +%s)
    local elapsed=$(( end_time - start_time ))
    echo "  FAILED after ${elapsed}s"
    write_marker "04_${grid_km}km_failed"
    echo "  [WARNING] Continuing despite 04 failure"
  fi

  # ── Stage 04b: 诊断 ─────────────────────────────────────────────
  run_stage "04b_${grid_km}km" "04b_recover_diagnostics.R" \
    "Recover diagnostics and species summaries"

  # ── Stage 05-07: 后处理 ─────────────────────────────────────────
  run_stage "05_${grid_km}km" "05_postprocess_diversity.R" \
    "Post-process diversity metrics + driver regression"

  run_stage "06_${grid_km}km" "06_figures_publication.R" \
    "Publication figures"

  run_stage "06b_${grid_km}km" "06b_regenerate_driver_plots.R" \
    "Driver + RF importance plots"

  run_stage "09_${grid_km}km" "09_extended_analyses.R" \
    "Extended analyses (varpart + RF)"

  run_stage "12_${grid_km}km" "12_homogenization_spatiotemporal_maps.R" \
    "Homogenization spatio-temporal maps"

  run_stage "14_${grid_km}km" "14_species_trait_regression.R" \
    "Species-trait regression (Q4)"

  echo ""
  echo "############################################################"
  echo "[$(timestamp)] Pipeline COMPLETE for ${grid_km}km grid"
  echo "############################################################"
}

# ══════════════════════════════════════════════════════════════════════
# 主流程：先 100km，再 10km
# ══════════════════════════════════════════════════════════════════════

echo "============================================================"
echo "[$(timestamp)] v3 STABLE pipeline — 100km then 10km"
echo "  Project:   $PROJECT_ROOT"
echo "  Log:       $PIPELINE_LOG"
echo "  Host:      $(hostname)"
echo "  tmux:      ${TMUX:-not in tmux}"
echo "============================================================"

write_marker "pipeline_start"

# ── Phase 1: 100km grid（可复用 v2 数据）────────────────────────────
run_pipeline_for_grid 100 "v3_full_200sp_ar1_spatial" "0"

# ── Phase 2: 10km grid（需要重新提取环境变量）────────────────────────
run_pipeline_for_grid 10 "v3_full_200sp_ar1_spatial" "0"

# ── Phase 3: 敏感性分析 ────────────────────────────────────────────
run_stage "15" "15_sensitivity_3yr_window.R" \
  "3-year window sensitivity analysis"

run_stage "15b" "15b_sensitivity_breeding_season.R" \
  "Breeding season sensitivity analysis"

write_marker "pipeline_all_done"
