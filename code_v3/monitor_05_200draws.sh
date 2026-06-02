#!/usr/bin/env bash
set -uo pipefail

# 监控 05_200draws (PID 914407)
PID=2108252
PROJECT_ROOT=~/bird_dynamic_occupancy_analysis
SCRIPT_DIR=$PROJECT_ROOT/code_v3
LOG_DIR=$PROJECT_ROOT/logs_v3
CHECK_INTERVAL=120

echo "[$(date)] 监控 PID $PID (05_postprocess_diversity_200draws)"
echo "[$(date)] 完成后自动运行 06-13 (200 draws pipeline)"

while kill -0 $PID 2>/dev/null; do
  sleep $CHECK_INTERVAL
done

echo "[$(date)] PID $PID 已结束"

# 检查关键输出
RUN_LABEL="v3_full_200sp_ar1_spatial_200draws"
if [ -f "$PROJECT_ROOT/data/derived_v3/psi_samples_thinned_${RUN_LABEL}.rds" ]; then
  echo "[$(date)] psi_samples_thinned_200draws 已生成"
else
  echo "[$(date)] 警告: psi_samples_thinned_200draws 未找到"
fi

echo "[$(date)] 启动 200 draws 06-13 pipeline"
cd $PROJECT_ROOT
bash $SCRIPT_DIR/run_pipeline_200draws.sh 2>&1 | tee $LOG_DIR/pipeline_200draws_auto_$(date +%Y%m%d_%H%M%S).log

echo "[$(date)] 200 draws pipeline 完成"
