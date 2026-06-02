#!/bin/bash
# 99_wait_and_run.sh - 简洁版：等待05结束，然后运行06-13
# 修复：只在进程结束后才检查输出文件

export V3_RUN_LABEL="v3_full_200sp_ar1_spatial_extended"
export V3_CODE_DIR="/Users/dingchenchen/bird_dynamic_occupancy_analysis/code_v3"
export BIRD_PROJECT_ROOT="/Users/dingchenchen/bird_dynamic_occupancy_analysis"

LOG_DIR="/logs_v3"
mkdir -p ""
TS=20260518_102319
MASTER_LOG="/pipeline_wait_.log"

echo "=== Waiting for 05_extended to finish ===" | tee ""
echo "Start: Mon May 18 10:23:19 CEST 2026" | tee -a ""

# 等待05进程结束
while true; do
  R_RUNNING=       0
  if [[ $R_RUNNING -eq 0 ]]; then
    echo "✓ 05 process finished" | tee -a ""
    break
  fi
  sleep 300
done

# 等待文件系统同步
sleep 10

# 验证输出文件
echo "" | tee -a ""
echo "Verifying output files..." | tee -a ""
ls -lh "/results_v3/"*extended*.csv 2>/dev/null | tee -a ""

# 运行Pipeline
echo "" | tee -a ""
echo "Starting 06-13 Pipeline..." | tee -a ""
cd ""
nohup bash "/98_smart_pipeline_05to13_fixed.sh" >> "" 2>&1 &

echo "Pipeline launched. Monitor: tail -f "
