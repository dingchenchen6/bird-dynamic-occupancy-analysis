#!/usr/bin/env bash
## 92_wait_upload_and_launch_server.sh
##
## 本地守护脚本：等待大文件上传结束，然后在服务器端自动
## 1. 解压 bundle
## 2. 设置 AVONET 路径
## 3. 使用 nohup 启动 v3 主流程

set -euo pipefail

UPLOAD_PID="${1:-0}"
SSH_PASSWORD="${2:?Need SSH password}"
REMOTE_HOST="${3:-dingchenchen@162.105.149.23}"
TASK_ROOT="${4:-/home/dingchenchen/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3}"
RUN_LABEL="${5:-v3_full_200sp_ar1_spatial}"
LOCAL_BUNDLE_PATH="${6:-/tmp/bird_dynamic_occupancy_analysis_v3_bundle.tar}"
REMOTE_BUNDLE="/home/dingchenchen/tmp/bird_dynamic_occupancy_analysis_v3_bundle.tar"
REMOTE_AVONET="${TASK_ROOT}/data/external/AVONET1_BirdLife.csv"

if [[ ! -f "${LOCAL_BUNDLE_PATH}" ]]; then
  echo "[watcher] ERROR: local bundle not found: ${LOCAL_BUNDLE_PATH}" >&2
  exit 1
fi

LOCAL_BUNDLE_SIZE="$(stat -f '%z' "${LOCAL_BUNDLE_PATH}")"
echo "[watcher] Monitoring remote bundle transfer at $(date '+%Y-%m-%d %H:%M:%S')"
echo "[watcher] Local bundle size: ${LOCAL_BUNDLE_SIZE}"

export SSHPASS="${SSH_PASSWORD}"

while true; do
  REMOTE_STATUS="$(
    /opt/homebrew/bin/sshpass -e ssh -o StrictHostKeyChecking=no "${REMOTE_HOST}" \
      "BUNDLE='${REMOTE_BUNDLE}'; AVONET='${REMOTE_AVONET}'; \
       if [ -f \"\$BUNDLE\" ]; then stat -c '%s' \"\$BUNDLE\" 2>/dev/null || stat -f '%z' \"\$BUNDLE\"; else echo 0; fi; \
       if [ -f \"\$AVONET\" ]; then echo avonet_ready; else echo avonet_pending; fi" \
      2>/dev/null || true
  )"

  REMOTE_BUNDLE_SIZE="$(printf '%s\n' "${REMOTE_STATUS}" | sed -n '1p' | tr -dc '0-9')"
  AVONET_FLAG="$(printf '%s\n' "${REMOTE_STATUS}" | sed -n '2p')"
  REMOTE_BUNDLE_SIZE="${REMOTE_BUNDLE_SIZE:-0}"

  echo "[watcher] Remote bundle size: ${REMOTE_BUNDLE_SIZE} | ${AVONET_FLAG} | $(date '+%Y-%m-%d %H:%M:%S')"

  if [[ "${REMOTE_BUNDLE_SIZE}" == "${LOCAL_BUNDLE_SIZE}" && "${AVONET_FLAG}" == "avonet_ready" ]]; then
    break
  fi

  sleep 60
done

echo "[watcher] Remote upload complete. Launching remote unpack/run at $(date '+%Y-%m-%d %H:%M:%S')"

/opt/homebrew/bin/sshpass -e ssh -o StrictHostKeyChecking=no "${REMOTE_HOST}" bash <<EOF
set -euo pipefail
TASK_ROOT="${TASK_ROOT}"
RUN_LABEL="${RUN_LABEL}"
REMOTE_BUNDLE="${REMOTE_BUNDLE}"
REMOTE_AVONET="${REMOTE_AVONET}"

mkdir -p "\${TASK_ROOT}" "\${TASK_ROOT}/logs_v3" "\${TASK_ROOT}/data/external"
tar -xf "\${REMOTE_BUNDLE}" -C "\${TASK_ROOT}"
chmod +x "\${TASK_ROOT}/code_v3/90_server_run_v3_pipeline.sh" \
         "\${TASK_ROOT}/code_v3/91_local_sync_to_server.sh" \
         "\${TASK_ROOT}/code_v3/92_wait_upload_and_launch_server.sh"

RUN_LOG="\${TASK_ROOT}/logs_v3/server_run_\${RUN_LABEL}_\$(date +%Y%m%d_%H%M%S).log"
nohup env \
  BIRD_PROJECT_ROOT="\${TASK_ROOT}" \
  V3_CODE_DIR="\${TASK_ROOT}/code_v3" \
  V3_AVONET_PATH="\${REMOTE_AVONET}" \
  V3_RUN_LABEL="\${RUN_LABEL}" \
  bash "\${TASK_ROOT}/code_v3/90_server_run_v3_pipeline.sh" "\${TASK_ROOT}" "\${RUN_LABEL}" \
  > "\${RUN_LOG}" 2>&1 < /dev/null &

echo \$! > "\${TASK_ROOT}/logs_v3/current_server_pid.txt"
echo "\${RUN_LOG}" > "\${TASK_ROOT}/logs_v3/current_server_log.txt"
echo "STARTED_PID=\$(cat "\${TASK_ROOT}/logs_v3/current_server_pid.txt")"
echo "RUN_LOG=\$(cat "\${TASK_ROOT}/logs_v3/current_server_log.txt")"
EOF

echo "[watcher] Remote launch submitted at $(date '+%Y-%m-%d %H:%M:%S')"
