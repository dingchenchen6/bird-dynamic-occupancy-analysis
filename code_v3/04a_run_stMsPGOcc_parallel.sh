#!/usr/bin/env bash
## 04a_run_stMsPGOcc_parallel.sh  —  4链并行运行 stMsPGOcc（线程安全版）
##
## 线程安全策略:
##   1. 每链绑定独立 CPU 亲和组（taskset），避免 OMP 线程跨链 false sharing
##   2. 每链 R 进程独立内存空间，不共享任何文件写入
##   3. OMP 线程数 = (总核数 - 留4核给系统) / 链数，避免超额订阅
##   4. OMP_NUM_THREADS 与 taskset 绑定范围一致
##   5. 保存前 gc() + fsync，确保数据落盘
##   6. 链完成检测加文件大小验证（≥10MB 才算有效）
##
## 用法：bash code_v3/04a_run_stMsPGOcc_parallel.sh

set -euo pipefail

PROJECT_ROOT="${BIRD_PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
SCRIPT_DIR="$PROJECT_ROOT/code_v3"
LOG_DIR="$PROJECT_ROOT/logs_v3"
mkdir -p "$LOG_DIR"

N_CHAINS="${V3_N_CHAINS:-4}"

# ── 线程安全: 动态计算 OMP 线程数 ──────────────────────────────────
# 留 4 核给系统/IO，剩余平分给各链
TOTAL_CORES=$(nproc 2>/dev/null || echo 256)
RESERVED=4
N_OMP_PER_CHAIN="${V3_N_OMP_PER_CHAIN:-$(( (TOTAL_CORES - RESERVED) / N_CHAINS ))}"
# 确保不低于 8 且不高于 63
N_OMP_PER_CHAIN=$(( N_OMP_PER_CHAIN < 8 ? 8 : N_OMP_PER_CHAIN ))
N_OMP_PER_CHAIN=$(( N_OMP_PER_CHAIN > 63 ? 63 : N_OMP_PER_CHAIN ))

echo "═══════════════════════════════════════════════════════════"
echo " stMsPGOcc 4-Chain Parallel (Thread-Safe)"
echo " Project:       $PROJECT_ROOT"
echo " Total cores:   $TOTAL_CORES"
echo " Chains:        $N_CHAINS"
echo " OMP/chain:     $N_OMP_PER_CHAIN"
echo " Total OMP:     $((N_CHAINS * N_OMP_PER_CHAIN))"
echo " CPU affinity:  taskset binding per chain"
echo " Start:         $(date '+%Y-%m-%d %H:%M:%S')"
echo "═══════════════════════════════════════════════════════════"

# ── 断点续跑: 检查已完成的链 ──────────────────────────────────────
SKIP_CHAINS=()
for CHAIN_ID in $(seq 1 "$N_CHAINS"); do
  CHAIN_FILE=$(find "$PROJECT_ROOT/data/derived_v3" -name "stMsPGOcc_fit_*_chain${CHAIN_ID}.rds" 2>/dev/null | head -1)
  if [[ -n "$CHAIN_FILE" ]]; then
    FILE_SIZE=$(stat -c%s "$CHAIN_FILE" 2>/dev/null || stat -f%z "$CHAIN_FILE" 2>/dev/null || echo "0")
    if [[ "$FILE_SIZE" -gt 10485760 ]]; then
      echo "  Chain $CHAIN_ID already exists (${FILE_SIZE} bytes) — skipping"
      SKIP_CHAINS+=("$CHAIN_ID")
    else
      echo "  Chain $CHAIN_ID file too small (${FILE_SIZE} bytes) — will re-run"
      rm -f "$CHAIN_FILE"
    fi
  fi
done

# ── 启动并行链（带 CPU 亲和绑定）─────────────────────────────────
PIDS=()
CHAIN_IDS=()

for CHAIN_ID in $(seq 1 "$N_CHAINS"); do
  if [[ " ${SKIP_CHAINS[*]} " == *" $CHAIN_ID "* ]]; then
    continue
  fi

  # CPU 亲和: 每链绑定到独立的核心组
  # 链1: 核心 0-62, 链2: 核心 63-125, 链3: 核心 126-188, 链4: 核心 189-251
  CPU_START=$(( (CHAIN_ID - 1) * N_OMP_PER_CHAIN ))
  CPU_END=$(( CPU_START + N_OMP_PER_CHAIN - 1 ))
  # 不超过最大核心数
  CPU_END=$(( CPU_END >= TOTAL_CORES ? TOTAL_CORES - 1 : CPU_END ))
  AFFINITY="${CPU_START}-${CPU_END}"

  CHAIN_LOG="${LOG_DIR}/04_chain${CHAIN_ID}_$(date +%Y%m%d_%H%M%S).log"
  echo ""
  echo "  [Chain $CHAIN_ID] CPU affinity: $AFFINITY, OMP=$N_OMP_PER_CHAIN"
  echo "  [Chain $CHAIN_ID] Log: $CHAIN_LOG"

  # 使用 taskset 绑定 CPU + 设置 OMP_NUM_THREADS
  V3_CHAIN_ID="$CHAIN_ID" \
  V3_N_OMP_THREADS="$N_OMP_PER_CHAIN" \
  OMP_NUM_THREADS="$N_OMP_PER_CHAIN" \
  OMP_PROC_BIND=close \
  OMP_PLACES=threads \
  taskset -c "$AFFINITY" \
    Rscript "${SCRIPT_DIR}/04_run_stMsPGOcc_main.R" >"$CHAIN_LOG" 2>&1 &

  PIDS+=($!)
  CHAIN_IDS+=("$CHAIN_ID")
done

if [[ ${#PIDS[@]} -eq 0 ]]; then
  echo ""
  echo "  All chains already completed. Nothing to run."
  echo "═══════════════════════════════════════════════════════════"
  exit 0
fi

# ── 等待所有链完成 ────────────────────────────────────────────────
echo ""
echo "  Waiting for ${#PIDS[@]} chain(s)..."
echo ""

FAILED=0
for i in "${!PIDS[@]}"; do
  PID="${PIDS[$i]}"
  CHAIN="${CHAIN_IDS[$i]}"
  if wait "$PID"; then
    echo "  [Chain $CHAIN] Done (PID $PID)"
    # 验证输出文件
    CHAIN_FILE=$(find "$PROJECT_ROOT/data/derived_v3" -name "stMsPGOcc_fit_*_chain${CHAIN}.rds" 2>/dev/null | head -1)
    if [[ -n "$CHAIN_FILE" ]]; then
      FSIZE=$(stat -c%s "$CHAIN_FILE" 2>/dev/null || stat -f%z "$CHAIN_FILE" 2>/dev/null || echo "0")
      if [[ "$FSIZE" -gt 10485760 ]]; then
        echo "  [Chain $CHAIN] Output verified: $(basename $CHAIN_FILE) ($(($FSIZE / 1048576)) MB)"
      else
        echo "  [Chain $CHAIN] WARNING: Output file too small ($FSIZE bytes)"
        FAILED=$((FAILED + 1))
      fi
    else
      echo "  [Chain $CHAIN] WARNING: No output file found"
      FAILED=$((FAILED + 1))
    fi
  else
    echo "  [Chain $CHAIN] FAILED (PID $PID, exit $?)"
    FAILED=$((FAILED + 1))
  fi
done

echo ""
echo "═══════════════════════════════════════════════════════════"
if [[ $FAILED -eq 0 ]]; then
  echo " All $N_CHAINS chains completed successfully"
  echo " Next: Run 04b_recover_diagnostics.R to combine chains"
else
  echo " WARNING: $FAILED chain(s) failed!"
  echo " Check logs: $LOG_DIR/04_chain*.log"
fi
echo " Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo "═══════════════════════════════════════════════════════════"

exit $FAILED
