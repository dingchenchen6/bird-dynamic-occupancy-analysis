#!/bin/bash
# 92_monitor_stage04.sh - 监控Stage 04运行状态
PROJ=~/bird_dynamic_occupancy_analysis
LOG_DIR="$PROJ/logs_v3"

echo "=========================================="
echo "Stage 04 Monitor - $(date)"
echo "=========================================="

# 1. 进程状态
echo ""
echo "【进程状态】"
ps aux | grep "04_run_stMsPGOcc_main.R" | grep -v grep | awk "{printf \"Chain PID=%s CPU=%s%% RSS=%.0fMB Time=%s\n\", \$2, \$3, \$6/1024, \$10}"

# 2. 日志进度
echo ""
echo "【日志进度】"
for i in 1 2 3 4; do
  LOG="$LOG_DIR/04_chain${i}_"*.log
  LATEST=$(ls -t $LOG 2>/dev/null | head -1)
  if [ -f "$LATEST" ]; then
    LINES=$(wc -l < "$LATEST")
    SIZE=$(stat -c%s "$LATEST" 2>/dev/null || stat -f%z "$LATEST")
    LAST_MOD=$(stat -c%Y "$LATEST" 2>/dev/null || stat -f%m "$LATEST")
    NOW=$(date +%s)
    AGE_MIN=$(( (NOW - LAST_MOD) / 60 ))
    echo "Chain $i: ${LINES} lines, ${SIZE} bytes, updated ${AGE_MIN} min ago"
    # 检查是否有segfault
    if grep -q "segfault\|Segmentation\|memory not mapped" "$LATEST" 2>/dev/null; then
      echo "  ⚠️  SEGFAULT detected!"
    fi
  else
    echo "Chain $i: no log"
  fi
done

# 3. 结果文件
echo ""
echo "【结果文件】"
ls -la "$PROJ/results_v3/"*.rds 2>/dev/null || echo "No RDS files yet"

# 4. 系统资源
echo ""
echo "【系统资源】"
free -h | grep Mem | awk "{printf \"Memory: used=%s free=%s\n\", \$3, \$4}"
df -h "$PROJ" | tail -1 | awk "{printf \"Disk: %s used, %s free\n\", \$3, \$4}"
