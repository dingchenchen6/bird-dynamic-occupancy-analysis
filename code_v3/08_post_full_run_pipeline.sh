#!/usr/bin/env bash
## 08_post_full_run_pipeline.sh  —  v3 全量运行后流水线
##
## 在 stMsPGOcc 全量拟合完成后，依次运行 stage 5-7 脚本。
## 用法：bash code_v3/08_post_full_run_pipeline.sh [RUN_LABEL]
##
## 阶段：
##   Stage 5:  后处理多样性 + 性状扩展 + 驱动因子
##   Stage 5b: RF 排列重要性
##   Stage 6:  出版级图表
##   Stage 6b: 驱动因子 + RF 图
##   Stage 7:  稿件 + 文档

set -euo pipefail

# ── 配置 ───────────────────────────────────────────────────────────────
RUN_LABEL="${1:-v3_full_200sp_ar1_spatial}"
export V3_RUN_LABEL="$RUN_LABEL"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$(dirname "$SCRIPT_DIR")/logs_v3"
mkdir -p "$LOG_DIR"

# R 可执行文件
RSCRIPT="$(which Rscript 2>/dev/null || echo /usr/local/bin/Rscript)"

echo "═══════════════════════════════════════════════════"
echo " v3 post-full-run pipeline"
echo " Run label: $RUN_LABEL"
echo " Started:   $(date '+%Y-%m-%d %H:%M:%S')"
echo "═══════════════════════════════════════════════════"

run_stage() {
  local stage_num="$1"
  local script="$2"
  local desc="$3"
  local log_file="${LOG_DIR}/${script%.R}_${RUN_LABEL}_$(date +%Y%m%d_%H%M%S).log"

  echo ""
  echo "── Stage ${stage_num}: ${desc} ──"
  echo "   Script: ${script}"
  echo "   Log:    ${log_file}"

  local start_time=$(date +%s)
  if "$RSCRIPT" "${SCRIPT_DIR}/${script}" 2>&1 | tee "$log_file"; then
    local end_time=$(date +%s)
    local elapsed=$(( end_time - start_time ))
    echo "   ✅ Done in ${elapsed}s"
  else
    local end_time=$(date +%s)
    local elapsed=$(( end_time - start_time ))
    echo "   ❌ FAILED after ${elapsed}s (exit code: $?)"
    echo "   Check log: ${log_file}"
    exit 1
  fi
}

# ── Stage 5: 性状扩展 ─────────────────────────────────────────────────
run_stage "5a" "03b_extend_traits.R" \
  "Extend traits (diet_specialization + habitat_breadth)"

# ── Stage 5: 后处理多样性 ──────────────────────────────────────────────
run_stage "5b" "05_postprocess_diversity.R" \
  "Post-process diversity metrics + brms driver regression"

# ── Stage 5c: Q4 性状回归 ──────────────────────────────────────────────
run_stage "5c" "14_species_trait_regression.R" \
  "Species-trait regression (Q4) with phylogenetic random effect"

# ── Stage 5d: 扩展分析（含 RF 重要性）──────────────────────────────────
run_stage "5d" "09_extended_analyses.R" \
  "Extended analyses (varpart + RF importance + Mann-Kendall)"

# ── Stage 5e: 同质化时空地图 ──────────────────────────────────────────
run_stage "5e" "12_homogenization_spatiotemporal_maps.R" \
  "Homogenization spatio-temporal maps"

# ── Stage 6: 出版级图表 ───────────────────────────────────────────────
run_stage "6a" "06_figures_publication.R" \
  "Publication figures (Nature style)"

run_stage "6b" "06b_regenerate_driver_plots.R" \
  "Driver + RF importance plots"

# ── Stage 7: 文档 ─────────────────────────────────────────────────────
run_stage "7a" "07_render_manuscript.R" \
  "Render workflow + manuscript markdown"

echo ""
echo "═══════════════════════════════════════════════════"
echo " v3 post-full-run pipeline COMPLETE"
echo " Run label: $RUN_LABEL"
echo " Finished:  $(date '+%Y-%m-%d %H:%M:%S')"
echo "═══════════════════════════════════════════════════"
