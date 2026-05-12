#!/usr/bin/env bash
## 90_server_run_v3_pipeline.sh  —  在服务器端运行 v3 动态占域分析
##
## 用途 / Purpose:
## 1. 统一在服务器端执行 code_v3 的前处理、主模型与后处理
## 2. stMsPGOcc 模型使用 4链并行（04a），加速 ~6-8x
## 3. 每个阶段独立写日志，便于排错与追踪进度
##
## 典型用法 / Example:
##   bash code_v3/90_server_run_v3_pipeline.sh \
##     /home/dingchenchen/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3 \
##     v3_full_200sp_ar1_spatial
##
## 并行设置:
##   V3_PARALLEL=1 (default): 4链并行, 每链32 OMP线程
##   V3_PARALLEL=0: 4链串行, 兼容模式

set -euo pipefail

PROJECT_ROOT="${1:-/home/dingchenchen/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3}"
RUN_LABEL="${2:-v3_full_200sp_ar1_spatial}"
PILOT_FLAG="${V3_PILOT:-0}"
PARALLEL_FLAG="${V3_PARALLEL:-1}"

export BIRD_PROJECT_ROOT="$PROJECT_ROOT"
export V3_CODE_DIR="$PROJECT_ROOT/code_v3"
export V3_RUN_LABEL="$RUN_LABEL"
export V3_PILOT="$PILOT_FLAG"

SCRIPT_DIR="$V3_CODE_DIR"
LOG_DIR="$PROJECT_ROOT/logs_v3"
mkdir -p "$LOG_DIR" "$PROJECT_ROOT/data/derived_v3" "$PROJECT_ROOT/results_v3" "$PROJECT_ROOT/figures_v3"

RSCRIPT="$(command -v Rscript || true)"
if [[ -z "$RSCRIPT" ]]; then
  echo "[server-run] ERROR: Rscript not found in PATH" >&2
  exit 1
fi

timestamp() {
  date '+%Y-%m-%d %H:%M:%S'
}

run_stage() {
  local stage_id="$1"
  local script_name="$2"
  local desc="$3"
  local log_file="${LOG_DIR}/${stage_id}_${script_name%.R}_${RUN_LABEL}_$(date +%Y%m%d_%H%M%S).log"

  echo ""
  echo "============================================================"
  echo "[server-run] Stage ${stage_id} :: ${desc}"
  echo "[server-run] Script: ${script_name}"
  echo "[server-run] Log:    ${log_file}"
  echo "[server-run] Start:  $(timestamp)"
  echo "============================================================"

  if "$RSCRIPT" "${SCRIPT_DIR}/${script_name}" >"$log_file" 2>&1; then
    echo "[server-run] Stage ${stage_id} completed: $(timestamp)"
  else
    echo "[server-run] Stage ${stage_id} failed: $(timestamp)" >&2
    echo "[server-run] Check log: ${log_file}" >&2
    exit 1
  fi
}

echo "============================================================"
echo "[server-run] v3 dynamic occupancy pipeline"
echo "[server-run] Project root: ${PROJECT_ROOT}"
echo "[server-run] Run label:    ${RUN_LABEL}"
echo "[server-run] Pilot mode:   ${PILOT_FLAG}"
echo "[server-run] Parallel:     ${PARALLEL_FLAG}"
echo "[server-run] Started:      $(timestamp)"
echo "============================================================"

## Stage 01-03: 前处理（串行）
run_stage "01" "01_merge_birdwatch_ebird.R" "Merge and re-deduplicate observation records"
run_stage "02" "02_build_survey_history.R" "Build breeding-season survey history"
run_stage "03" "03_prepare_environment.R" "Prepare baseline environment and detection covariates"
run_stage "03b" "03b_extend_traits.R" "Extend functional trait table"
run_stage "03c" "03c_prepare_climate_change.R" "Prepare climate-change covariates"
run_stage "03d" "03d_prepare_landuse_change.R" "Prepare land-use-change covariates"
run_stage "03e" "03e_prepare_hfi_change.R" "Prepare human footprint change covariates"

## Stage 04: stMsPGOcc 模型拟合
if [[ "$PARALLEL_FLAG" == "1" ]]; then
  echo ""
  echo "============================================================"
  echo "[server-run] Stage 04 :: stMsPGOcc 4-chain PARALLEL fit"
  echo "[server-run] OMP threads/chain: ${V3_N_OMP_PER_CHAIN:-32}"
  echo "[server-run] Start:             $(timestamp)"
  echo "============================================================"

  export V3_N_OMP_PER_CHAIN="${V3_N_OMP_PER_CHAIN:-32}"
  export V3_N_CHAINS="${V3_N_CHAINS:-4}"

  if bash "${SCRIPT_DIR}/04a_run_stMsPGOcc_parallel.sh" >"${LOG_DIR}/04_parallel_${RUN_LABEL}_$(date +%Y%m%d_%H%M%S).log" 2>&1; then
    echo "[server-run] Stage 04 parallel completed: $(timestamp)"
  else
    echo "[server-run] Stage 04 parallel FAILED: $(timestamp)" >&2
    exit 1
  fi
else
  run_stage "04" "04_run_stMsPGOcc_main.R" "Fit multi-species spatiotemporal occupancy model (sequential)"
fi

## Stage 04b: 诊断与链合并
run_stage "04b" "04b_recover_diagnostics.R" "Combine chains, recover diagnostics and species-level summaries"

## Stage 05-07: 后处理、图件、文稿
bash "${SCRIPT_DIR}/08_post_full_run_pipeline.sh" "${RUN_LABEL}" >"${LOG_DIR}/08_post_full_run_pipeline_${RUN_LABEL}_$(date +%Y%m%d_%H%M%S).log" 2>&1

echo "============================================================"
echo "[server-run] v3 pipeline finished successfully"
echo "[server-run] Finished: $(timestamp)"
echo "============================================================"
