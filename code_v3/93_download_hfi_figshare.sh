#!/usr/bin/env bash
## 93_download_hfi_figshare.sh  —  从Figshare下载HFI 2000-2018数据
##
## 数据集: Mu et al. 2022, An annual global terrestrial Human Footprint dataset from 2000 to 2018
## DOI: 10.6084/m9.figshare.16571064
##
## 用法:
##   bash code_v3/93_download_hfi_figshare.sh [下载目录]

set -euo pipefail

DOWNLOAD_DIR="${1:-/home/dingchenchen/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3/data/hfi_figshare}"
mkdir -p "$DOWNLOAD_DIR"

echo "========================================"
echo "HFI Figshare 下载脚本"
echo "目标目录: $DOWNLOAD_DIR"
echo "========================================"

# Figshare base URL (通过DOI解析)
# 文件命名格式: hfp2000.zip - hfp2018.zip
BASE_URL="https://figshare.com/ndownloader/files"

# 根据Figshare API，文件ID列表（需要确认）
# 备用：直接通过文章页面获取下载链接

# 方法1: 直接下载（如果URL模式已知）
# 方法2: 通过wget递归下载页面中的zip文件

# 先尝试通过DOI landing page获取下载链接
echo "[1/3] 获取Figshare下载页面..."

# Figshare文章页面
FIGSHARE_ARTICLE="https://doi.org/10.6084/m9.figshare.16571064"

# 使用curl跟随重定向获取实际下载链接
echo "[2/3] 解析下载链接..."

# 创建下载函数
download_year() {
  local year=$1
  local url=$2
  local output="$DOWNLOAD_DIR/hfp${year}.zip"
  
  if [[ -f "$output" ]]; then
    echo "  [跳过] hfp${year}.zip 已存在"
    return 0
  fi
  
  echo "  [下载] hfp${year}.zip ..."
  if curl -L --max-time 300 -o "$output" "$url" 2>/dev/null; then
    local size=$(du -h "$output" 2>/dev/null | cut -f1)
    echo "  [完成] hfp${year}.zip ($size)"
    return 0
  else
    echo "  [失败] hfp${year}.zip"
    rm -f "$output"
    return 1
  fi
}

# 方法：通过Figshare API获取文件列表
echo "[3/3] 通过Figshare API获取文件列表并下载..."

# Figshare article ID: 16571064
ARTICLE_ID="16571064"

# 获取文件列表
FILES_JSON=$(curl -s "https://api.figshare.com/v2/articles/${ARTICLE_ID}/files")

# 解析并下载每个文件
echo "$FILES_JSON" | python3 -c "
import sys, json
try:
    files = json.load(sys.stdin)
    for f in files:
        name = f.get('name', '')
        url = f.get('download_url', '')
        if name.endswith('.zip') and url:
            print(f'{name}|{url}')
except Exception as e:
    print(f'ERROR: {e}', file=sys.stderr)
" 2>/dev/null | while IFS='|' read -r name url; do
  [[ -z "$name" ]] && continue
  year=$(echo "$name" | grep -oP 'hfp\K\d{4}' || echo "")
  [[ -z "$year" ]] && continue
  
  output="$DOWNLOAD_DIR/$name"
  if [[ -f "$output" ]]; then
    echo "  [跳过] $name 已存在"
    continue
  fi
  
  echo "  [下载] $name ..."
  if curl -L --max-time 600 -o "$output" "$url" 2>/dev/null; then
    size=$(du -h "$output" 2>/dev/null | cut -f1)
    echo "  [完成] $name ($size)"
  else
    echo "  [失败] $name"
    rm -f "$output"
  fi
done

echo ""
echo "========================================"
echo "下载完成！"
echo "文件列表:"
ls -lh "$DOWNLOAD_DIR"/*.zip 2>/dev/null || echo "  无zip文件"
echo "========================================"

# 统计
TOTAL_FILES=$(ls "$DOWNLOAD_DIR"/*.zip 2>/dev/null | wc -l)
echo "总计: $TOTAL_FILES 个zip文件"
