#!/usr/bin/env bash
## 99_launch_pipeline_extended.sh
## 启动06-13 Pipeline for extended版本 (400 draws)
## 设置正确的V3_RUN_LABEL环境变量

export V3_RUN_LABEL="v3_full_200sp_ar1_spatial_extended"
export V3_CODE_DIR="$HOME/bird_dynamic_occupancy_analysis/code_v3"
export BIRD_PROJECT_ROOT="$HOME/bird_dynamic_occupancy_analysis"

echo "=========================================="
echo " v3 Pipeline Launcher for EXTENDED (400 draws)"
echo "=========================================="
echo " V3_RUN_LABEL = $V3_RUN_LABEL"
echo " Start time: $(date +%Y-%m-%d %H:%M:%S)"
echo "=========================================="

cd "$BIRD_PROJECT_ROOT" || exit 1

# 检查05_extended输出文件是否存在
OUTPUT_PATTERN="results_v3/*${V3_RUN_LABEL}*.csv"
if ls $OUTPUT_PATTERN 1>/dev/null 2>&1; then
    echo "✓ Found 05_extended output files:"
    ls -lh $OUTPUT_PATTERN | head -5
else
    echo "⚠️  Warning: No output files found matching pattern: $OUTPUT_PATTERN"
    echo "   Waiting 60s then checking again..."
    sleep 60
fi

# 启动Pipeline监控脚本
echo ""
echo "Starting smart pipeline..."
exec "$V3_CODE_DIR/98_smart_pipeline_05to13.sh"
