#!/usr/bin/env bash
set -uo pipefail

PID=3311611
PROJECT_ROOT=~/bird_dynamic_occupancy_analysis
SCRIPT_DIR=$PROJECT_ROOT/code_v3
LOG_DIR=$PROJECT_ROOT/logs_v3
CHECK_INTERVAL=120

echo "[$(date)] 监控 PID $PID (05_postprocess_diversity)"
echo "[$(date)] 完成后自动运行 06-13 pipeline"

while kill -0 $PID 2>/dev/null; do
  sleep $CHECK_INTERVAL
done

EXIT_CODE=$(ps -p $PID -o exit= 2>/dev/null | awk '{print $1}')
echo "[$(date)] PID $PID 已结束"

# 检查05是否生成了关键输出
if [[ -f "$PROJECT_ROOT/data/derived_v3/psi_samples_thinned_v3_full_200sp_ar1_spatial.rds" ]]; then
  echo "[$(date)] psi_samples_thinned 已生成"
else
  echo "[$(date)] 警告: psi_samples_thinned 未找到，Stage 12 可能失败"
fi

# 检查05日志是否报错
LOG05=$(ls -t $LOG_DIR/05_rerun_fixed2_*.log 2>/dev/null | head -1)
if [[ -n "$LOG05" ]]; then
  if grep -q 'Error' "$LOG05"; then
    echo "[$(date)] 05 日志中有 Error，但继续运行 06-13"
  else
    echo "[$(date)] 05 日志无 Error"
  fi
fi

echo "[$(date)] 启动 06-13 pipeline"
cd $PROJECT_ROOT
export V3_CODE_DIR=$SCRIPT_DIR
export BIRD_PROJECT_ROOT=$PROJECT_ROOT
bash $SCRIPT_DIR/98_smart_pipeline_05to13.sh 2>&1 | tee $LOG_DIR/pipeline_06to13_after05_$(date +%Y%m%d_%H%M%S).log

echo "[$(date)] 06-13 pipeline 完成"
