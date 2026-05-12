#!/usr/bin/env bash
## 91_local_sync_to_server.sh  —  从本地同步 v3 项目到服务器
##
## 说明 / Notes:
## - 此脚本在本地执行
## - 只同步 90_sync_manifest_v3.txt 中列出的必要目录与文件
## - 不直接启动分析；同步完成后可再调用 90_server_run_v3_pipeline.sh

set -euo pipefail

LOCAL_ROOT="${1:-$HOME/Documents/New project/bird_dynamic_occupancy_analysis}"
REMOTE_HOST="${2:-dingchenchen@162.105.149.23}"
REMOTE_ROOT="${3:-/home/dingchenchen/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3}"
MANIFEST="${LOCAL_ROOT}/code_v3/90_sync_manifest_v3.txt"

if [[ ! -f "$MANIFEST" ]]; then
  echo "[local-sync] ERROR: manifest not found: $MANIFEST" >&2
  exit 1
fi

mkdir -p "${LOCAL_ROOT}/logs_v3"

echo "[local-sync] Local root:  ${LOCAL_ROOT}"
echo "[local-sync] Remote host: ${REMOTE_HOST}"
echo "[local-sync] Remote root: ${REMOTE_ROOT}"

while IFS= read -r rel_path; do
  [[ -z "$rel_path" ]] && continue
  echo "[local-sync] Syncing: ${rel_path}"
  rsync -a "${LOCAL_ROOT}/${rel_path}" "${REMOTE_HOST}:${REMOTE_ROOT}/$(dirname "${rel_path}")/"
done < "$MANIFEST"

echo "[local-sync] Sync completed."
