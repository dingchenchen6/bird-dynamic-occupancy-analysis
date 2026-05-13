#!/usr/bin/env bash
## 98_smart_pipeline_05to13.sh — 智能连续运行管线
## 在 05 完成后自动运行 06-13，带输入检查

set -uo pipefail

# ── 配置 ──────────────────────────────────────────────
RUN_LABEL="${V3_RUN_LABEL:-v3_full_200sp_ar1_spatial}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
LOG_DIR="$PROJECT_ROOT/logs_v3"
mkdir -p "$LOG_DIR"

TS="$(date +%Y%m%d_%H%M%S)"
MASTER_LOG="$LOG_DIR/pipeline_smart_${RUN_LABEL}_${TS}.log"

# R 可执行文件
RSCRIPT="$(which Rscript 2>/dev/null || echo /opt/R/4.5.3/lib64/R/bin/Rscript)"

echo "═══════════════════════════════════════════════════════" | tee -a "$MASTER_LOG"
echo " v3 Smart Pipeline: 05 → 13 (auto-continue)" | tee -a "$MASTER_LOG"
echo " Run label:  $RUN_LABEL" | tee -a "$MASTER_LOG"
echo " Started:    $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "$MASTER_LOG"
echo " Master log: $MASTER_LOG" | tee -a "$MASTER_LOG"
echo "═══════════════════════════════════════════════════════" | tee -a "$MASTER_LOG"

# ── 工具函数 ──────────────────────────────────────────────

# 运行 R 脚本（带重试）
run_r_stage() {
  local stage_num="$1"
  local script="$2"
  local desc="$3"
  local script_path="${SCRIPT_DIR}/${script}"
  local log_file="${LOG_DIR}/${stage_num}_${script%.R}_${TS}.log"

  # 设置环境变量
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

  # 06 需要检查输入文件
  if [[ "$stage_num" == "06" || "$stage_num" == "06b" ]]; then
    local input_file="$PROJECT_ROOT/results_v3/table_community_metrics_with_cri_${RUN_LABEL}.csv"
    if [[ ! -f "$input_file" ]]; then
      echo "   ⚠️  Input file not found for ${stage_num}: $input_file" | tee -a "$MASTER_LOG"
      echo "   Skipping ${stage_num}" | tee -a "$MASTER_LOG"
      return 1
    fi
    local lines=$(wc -l < "$input_file" 2>/dev/null || echo 0)
    if [[ $lines -lt 10 ]]; then
      echo "   ⚠️  Input file appears empty for ${stage_num}: ${lines} lines" | tee -a "$MASTER_LOG"
      echo "   Skipping ${stage_num}" | tee -a "$MASTER_LOG"
      return 1
    fi
    echo "   Input OK: $input_file (${lines} lines)" | tee -a "$MASTER_LOG"
  fi

  local start_time=$(date +%s)
  local exec_cmd="$RSCRIPT"
  if [[ "$script" == *.sh ]]; then
    exec_cmd="bash"
  fi
  if $exec_cmd "$script_path" > "$log_file" 2>&1; then
    local end_time=$(date +%s)
    local elapsed=$(( end_time - start_time ))
    local elapsed_min=$(( elapsed / 60 ))
    echo "   ✅ Done in ${elapsed}s (${elapsed_min}m)" | tee -a "$MASTER_LOG"
    return 0
  else
    local exit_code=$?
    local end_time=$(date +%s)
    local elapsed=$(( end_time - start_time ))
    echo "   ❌ FAILED after ${elapsed}s (exit code: ${exit_code})" | tee -a "$MASTER_LOG"
    echo "   Check log: ${log_file}" | tee -a "$MASTER_LOG"
    return 1
  fi
}

# ── 等待 05 完成 ──────────────────────────────────────────────
echo "" | tee -a "$MASTER_LOG"
echo "🔍 Waiting for 05 to complete..." | tee -a "$MASTER_LOG"

while true; do
  # 检查 05 R 进程是否还在
  R_RUNNING=$(ps aux | grep "05_postprocess_diversity.R" | grep -v grep | wc -l)
  
  if [[ $R_RUNNING -eq 0 ]]; then
    echo "✅ 05 R process finished (not found in ps)" | tee -a "$MASTER_LOG"
    break
  fi
  
  # 检查输出文件是否生成
  OUTPUT_FILE="$PROJECT_ROOT/results_v3/table_diversity_summary_${RUN_LABEL}.csv"
  if [[ -f "$OUTPUT_FILE" ]]; then
    LINES=$(wc -l < "$OUTPUT_FILE" 2>/dev/null || echo 0)
    if [[ $LINES -gt 100 ]]; then
      echo "✅ 05 output detected: $OUTPUT_FILE (${LINES} lines)" | tee -a "$MASTER_LOG"
      break
    fi
  fi
  
  echo "   ... still running ($(date '+%H:%M:%S'))" | tee -a "$MASTER_LOG"
  sleep 300  # 每 5 分钟检查一次
done

# 给一点时间让文件完全写入
sleep 10

# ── 运行 06-13 ──────────────────────────────────────────────
FAILED_STAGES=()

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

# Stage 8: 全量后处理流水线
if ! run_r_stage "08" "08_post_full_run_pipeline.sh" "Full post-run pipeline"; then
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

# ── 汇总 ──────────────────────────────────────────────
echo "" | tee -a "$MASTER_LOG"
echo "═══════════════════════════════════════════════════════" | tee -a "$MASTER_LOG"

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
echo "═══════════════════════════════════════════════════════" | tee -a "$MASTER_LOG"
