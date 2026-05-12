#!/usr/bin/env bash
## 95_server_hfi_and_pipeline.sh  —  服务器端HFI下载+v3 Pipeline运行
##
## 用法:
##   export V3_GRID_SIZE_KM="10"
##   export V3_MIN_VISITS="1"
##   export V3_MIN_SPECIES_GRID="3"
##   bash code_v3/95_server_hfi_and_pipeline.sh

set -euo pipefail

# 配置
PROJECT_ROOT="/home/dingchenchen/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3"
HFI_DIR="$PROJECT_ROOT/data/hfi_wcs"
CODE_DIR="$PROJECT_ROOT/code_v3"
LOG_DIR="$PROJECT_ROOT/logs_v3"

mkdir -p "$HFI_DIR" "$LOG_DIR"

# 环境变量
export BIRD_PROJECT_ROOT="$PROJECT_ROOT"
export V3_CODE_DIR="$CODE_DIR"
export V3_HFI_DIR="$HFI_DIR"

# 网格配置
GRID_SIZE_KM="${V3_GRID_SIZE_KM:-100}"
GRID_TAG=""
if [[ "$GRID_SIZE_KM" != "100" ]]; then
  GRID_TAG="_${GRID_SIZE_KM}km"
fi

RUN_LABEL="v3_full_200sp_ar1_spatial${GRID_TAG}"
export V3_RUN_LABEL="$RUN_LABEL"

echo "========================================"
echo "服务器 HFI 下载 + v3 Pipeline"
echo "========================================"
echo "项目目录: $PROJECT_ROOT"
echo "网格大小: ${GRID_SIZE_KM}km"
echo "运行标签: $RUN_LABEL"
echo "========================================"

# 函数：记录日志
log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# ========================================
# Phase 1: 下载HFI数据
# ========================================
echo ""
echo "========================================"
echo "Phase 1: HFI数据下载"
echo "========================================"

# 检查是否已有HFI数据
HFI_COUNT=$(ls "$HFI_DIR"/hfi_*.tif 2>/dev/null | wc -l)
if [[ $HFI_COUNT -ge 20 ]]; then
  log "HFI数据已存在 ($HFI_COUNT 个文件)，跳过下载"
else
  log "开始下载WCS HFI数据..."
  bash "$CODE_DIR/93_download_hfi_wcs.sh" "$HFI_DIR" \
    > "$LOG_DIR/hfi_download_${RUN_LABEL}.log" 2>&1
  log "HFI下载完成"
fi

# ========================================
# Phase 2: 运行v3 Pipeline
# ========================================
echo ""
echo "========================================"
echo "Phase 2: v3 Pipeline"
echo "========================================"

# 设置HFI环境变量供R脚本使用
export V3_HFI_WCS_DIR="$HFI_DIR"

# 运行完整pipeline
bash "$CODE_DIR/90_server_run_v3_pipeline.sh" "$PROJECT_ROOT" "$RUN_LABEL" \
  > "$LOG_DIR/pipeline_${RUN_LABEL}.log" 2>&1

echo ""
echo "========================================"
echo "全部完成！"
echo "========================================"
echo "HFI数据: $HFI_DIR"
echo "日志文件: $LOG_DIR"
echo "结果文件: $PROJECT_ROOT/results_v3"
echo "图件文件: $PROJECT_ROOT/figures_v3"
echo "========================================"
