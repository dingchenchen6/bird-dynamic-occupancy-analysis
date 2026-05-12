#!/usr/bin/env bash
## 98_nohup_run_pipeline_05to13.sh  —  v3 后台连续运行 05→13
##
## 用法：
##   nohup bash code_v3/98_nohup_run_pipeline_05to13.sh > logs_v3/nohup_pipeline_05to13_$(date +%Y%m%d_%H%M%S).log 2>&1 &
##
## 或：
##   bash code_v3/98_nohup_run_pipeline_05to13.sh &
##
## 阶段：
##   Stage 5:  后处理多样性
##   Stage 6:  出版级图表
##   Stage 7:  稿件渲染
##   Stage 8:  全量后处理流水线（03b, 05, 14, 09, 12, 06, 06b, 07）
##   Stage 10: PPTX 输出
##   Stage 11: DOCX 期刊文档
##   Stage 12: 同质化地图
##   Stage 13: 全球提案 DOCX

set -euo pipefail

# ── 配置 ───────────────────────────────────────────────────────────────
RUN_LABEL="${V3_RUN_LABEL:-v3_full_200sp_ar1_spatial}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
LOG_DIR="$PROJECT_ROOT/logs_v3"
mkdir -p "$LOG_DIR"

# 时间戳
TS="$(date +%Y%m%d_%H%M%S)"
MASTER_LOG="$LOG_DIR/pipeline_05to13_${RUN_LABEL}_${TS}.log"

# R 可执行文件
RSCRIPT="$(which Rscript 2>/dev/null || echo /usr/local/bin/Rscript)"

echo "═══════════════════════════════════════════════════════════════" | tee -a "$MASTER_LOG"
echo " v3 Pipeline: 05 → 13 (后台连续运行)" | tee -a "$MASTER_LOG"
echo " Run label:  $RUN_LABEL" | tee -a "$MASTER_LOG"
echo " Started:    $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "$MASTER_LOG"
echo " Master log: $MASTER_LOG" | tee -a "$MASTER_LOG"
echo "═══════════════════════════════════════════════════════════════" | tee -a "$MASTER_LOG"

# ── 运行单个 R 脚本 ────────────────────────────────────────────────────
run_r_stage() {
  local stage_num="$1"
  local script="$2"
  local desc="$3"
  local script_path="${SCRIPT_DIR}/${script}"
  local log_file="${LOG_DIR}/${stage_num}_${script%.R}_${TS}.log"
  # 设置 R 脚本需要的环境变量
  export V3_CODE_DIR="$SCRIPT_DIR"
  export BIRD_PROJECT_ROOT="$PROJECT_ROOT"

  echo "" | tee -a "$MASTER_LOG"
  echo "── Stage ${stage_num}: ${desc} ──" | tee -a "$MASTER_LOG"
  echo "   Script: ${script_path}" | tee -a "$MASTER_LOG"
  echo "   Log:    ${log_file}" | tee -a "$MASTER_LOG"

  if [[ ! -f "$script_path" ]]; then
    echo "   ⚠️  Script not found, skipping" | tee -a "$MASTER_LOG"
    return 0
  fi

  local start_time=$(date +%s)
  if "$RSCRIPT" "$script_path" > "$log_file" 2>&1; then
    local end_time=$(date +%s)
    local elapsed=$(( end_time - start_time ))
    local elapsed_min=$(( elapsed / 60 ))
    echo "   ✅ Done in ${elapsed}s (${elapsed_min}m)" | tee -a "$MASTER_LOG"
  else
    local exit_code=$?
    local end_time=$(date +%s)
    local elapsed=$(( end_time - start_time ))
    echo "   ❌ FAILED after ${elapsed}s (exit code: ${exit_code})" | tee -a "$MASTER_LOG"
    echo "   Check log: ${log_file}" | tee -a "$MASTER_LOG"
    return 1
  fi
}

# ── 运行 shell 脚本 ────────────────────────────────────────────────────
run_sh_stage() {
  local stage_num="$1"
  local script="$2"
  local desc="$3"
  local script_path="${SCRIPT_DIR}/${script}"
  local log_file="${LOG_DIR}/${stage_num}_${script%.sh}_${TS}.log"
  # 设置 R 脚本需要的环境变量
  export V3_CODE_DIR="$SCRIPT_DIR"
  export BIRD_PROJECT_ROOT="$PROJECT_ROOT"

  echo "" | tee -a "$MASTER_LOG"
  echo "── Stage ${stage_num}: ${desc} ──" | tee -a "$MASTER_LOG"
  echo "   Script: ${script_path}" | tee -a "$MASTER_LOG"
  echo "   Log:    ${log_file}" | tee -a "$MASTER_LOG"

  if [[ ! -f "$script_path" ]]; then
    echo "   ⚠️  Script not found, skipping" | tee -a "$MASTER_LOG"
    return 0
  fi

  local start_time=$(date +%s)
  if bash "$script_path" "$RUN_LABEL" > "$log_file" 2>&1; then
    local end_time=$(date +%s)
    local elapsed=$(( end_time - start_time ))
    local elapsed_min=$(( elapsed / 60 ))
    echo "   ✅ Done in ${elapsed}s (${elapsed_min}m)" | tee -a "$MASTER_LOG"
  else
    local exit_code=$?
    local end_time=$(date +%s)
    local elapsed=$(( end_time - start_time ))
    echo "   ❌ FAILED after ${elapsed}s (exit code: ${exit_code})" | tee -a "$MASTER_LOG"
    echo "   Check log: ${log_file}" | tee -a "$MASTER_LOG"
    return 1
  fi
}

# ── Pipeline 执行 ───────────────────────────────────────────────────────
FAILED_STAGES=()

# Stage 5: 后处理多样性（核心，最耗时）
if ! run_r_stage "05" "05_postprocess_diversity.R" "Post-process diversity metrics"; then
  FAILED_STAGES+=("05")
fi

# Stage 6: 出版级图表
if ! run_r_stage "06" "06_figures_publication.R" "Publication figures"; then
  FAILED_STAGES+=("06")
fi

# Stage 6b: 驱动因子图
if ! run_r_stage "06b" "06b_regenerate_driver_plots.R" "Driver + RF plots"; then
  FAILED_STAGES+=("06b")
fi

# Stage 7: 稿件渲染
if ! run_r_stage "07" "07_render_manuscript.R" "Render manuscript"; then
  FAILED_STAGES+=("07")
fi

# Stage 8: 全量后处理流水线（08 内部会运行 03b, 14, 09, 12）
if ! run_sh_stage "08" "08_post_full_run_pipeline.sh" "Full post-run pipeline"; then
  FAILED_STAGES+=("08")
fi

# Stage 10: PPTX
if ! run_r_stage "10" "10_finalize_with_dashline_pptx.R" "PPTX output"; then
  FAILED_STAGES+=("10")
fi

# Stage 11: DOCX 期刊
if ! run_r_stage "11" "11_render_journal_docx.R" "Journal DOCX"; then
  FAILED_STAGES+=("11")
fi

# Stage 12: 同质化地图
if ! run_r_stage "12" "12_homogenization_spatiotemporal_maps.R" "Homogenization maps"; then
  FAILED_STAGES+=("12")
fi

# Stage 13: 全球提案
if ! run_r_stage "13" "13_render_global_proposal_docx.R" "Global proposal DOCX"; then
  FAILED_STAGES+=("13")
fi

# ── 汇总 ────────────────────────────────────────────────────────────────
echo "" | tee -a "$MASTER_LOG"
echo "═══════════════════════════════════════════════════════════════" | tee -a "$MASTER_LOG"

if [[ ${#FAILED_STAGES[@]} -eq 0 ]]; then
  echo " ✅ v3 Pipeline 05→13 ALL COMPLETE" | tee -a "$MASTER_LOG"
else
  echo " ⚠️  v3 Pipeline completed with failures:" | tee -a "$MASTER_LOG"
  for s in "${FAILED_STAGES[@]}"; do
    echo "    - Stage ${s} FAILED" | tee -a "$MASTER_LOG"
  done
fi

echo " Run label:  $RUN_LABEL" | tee -a "$MASTER_LOG"
echo " Finished:   $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "$MASTER_LOG"
echo " Master log: $MASTER_LOG" | tee -a "$MASTER_LOG"
echo "═══════════════════════════════════════════════════════════════" | tee -a "$MASTER_LOG"
