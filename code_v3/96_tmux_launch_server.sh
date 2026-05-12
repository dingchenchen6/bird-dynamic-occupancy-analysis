#!/usr/bin/env bash
## 96_tmux_launch_server.sh  —  在服务器通过tmux启动HFI下载和Pipeline
##
## 用法:
##   ssh admin@162.105.149.23
##   bash /path/to/96_tmux_launch_server.sh

set -euo pipefail

# 配置
PROJECT_ROOT="/home/dingchenchen/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3"
CODE_DIR="$PROJECT_ROOT/code_v3"
LOG_DIR="$PROJECT_ROOT/logs_v3"

# 环境变量
export V3_GRID_SIZE_KM="${V3_GRID_SIZE_KM:-100}"
export V3_MIN_VISITS="${V3_MIN_VISITS:-1}"
export V3_MIN_SPECIES_GRID="${V3_MIN_SPECIES_GRID:-3}"
export V3_N_OMP_THREADS="${V3_N_OMP_THREADS:-16}"

GRID_TAG=""
if [[ "$V3_GRID_SIZE_KM" != "100" ]]; then
  GRID_TAG="_${V3_GRID_SIZE_KM}km"
fi

RUN_LABEL="v3_full_200sp_ar1_spatial${GRID_TAG}"

echo "========================================"
echo "Tmux 会话启动脚本"
echo "========================================"
echo "网格大小: ${V3_GRID_SIZE_KM}km"
echo "运行标签: $RUN_LABEL"
echo "========================================"

# 创建tmux会话: hfi_download
if tmux has-session -t hfi_download 2>/dev/null; then
  echo "[hfi_download] 会话已存在，附加..."
else
  echo "[hfi_download] 创建新会话..."
  tmux new-session -d -s hfi_download
  tmux send-keys -t hfi_download "cd $PROJECT_ROOT" C-m
  tmux send-keys -t hfi_download "echo '开始HFI数据下载...'" C-m
  tmux send-keys -t hfi_download "bash $CODE_DIR/93_download_hfi_wcs.sh $PROJECT_ROOT/data/hfi_wcs 2>&1 | tee $LOG_DIR/hfi_download_${RUN_LABEL}.log" C-m
fi

# 创建tmux会话: bird_pipeline
if tmux has-session -t bird_pipeline 2>/dev/null; then
  echo "[bird_pipeline] 会话已存在，附加..."
else
  echo "[bird_pipeline] 创建新会话..."
  tmux new-session -d -s bird_pipeline
  tmux send-keys -t bird_pipeline "cd $PROJECT_ROOT" C-m
  tmux send-keys -t bird_pipeline "echo '等待HFI下载完成...'" C-m
  tmux send-keys -t bird_pipeline "while [[ \$(ls $PROJECT_ROOT/data/hfi_wcs/hfi_*.tif 2>/dev/null | wc -l) -lt 5 ]]; do sleep 60; done" C-m
  tmux send-keys -t bird_pipeline "echo 'HFI数据就绪，开始Pipeline...'" C-m
  tmux send-keys -t bird_pipeline "export V3_HFI_WCS_DIR=$PROJECT_ROOT/data/hfi_wcs" C-m
  tmux send-keys -t bird_pipeline "bash $CODE_DIR/95_server_hfi_and_pipeline.sh 2>&1 | tee $LOG_DIR/pipeline_${RUN_LABEL}.log" C-m
fi

echo ""
echo "========================================"
echo "Tmux 会话已创建:"
echo "  - hfi_download: HFI数据下载"
echo "  - bird_pipeline: v3 Pipeline运行"
echo ""
echo "常用命令:"
echo "  tmux ls                    # 列出会话"
echo "  tmux attach -t hfi_download    # 查看下载进度"
echo "  tmux attach -t bird_pipeline   # 查看Pipeline进度"
echo "  tmux detach                # 分离会话 (Ctrl+B, D)"
echo "========================================"
