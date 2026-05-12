#!/usr/bin/env bash
## 93_download_hfi_wcs.sh  —  从WCS下载HFI 2001-2020数据
##
## 数据源: Wildlife Conservation Society (WCS)
## 网址: https://www.wcshumanfootprint.org/data-access
##
## 说明:
## - 2001-2014: HF3-1方法
## - 2015-2020: HF3-2方法 (推荐，有OSM数据)
## - 2000年数据缺失，可用2001年代替或联系作者
##
## 用法:
##   bash code_v3/93_download_hfi_wcs.sh [下载目录]

set -euo pipefail

DOWNLOAD_DIR="${1:-/home/dingchenchen/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3/data/hfi_wcs}"
mkdir -p "$DOWNLOAD_DIR"

echo "========================================"
echo "HFI WCS 下载脚本"
echo "数据源: Wildlife Conservation Society"
echo "目标目录: $DOWNLOAD_DIR"
echo "========================================"

# 下载函数
download_hfi() {
  local year=$1
  local url=$2
  local output="$DOWNLOAD_DIR/hfi_${year}.tif"
  
  if [[ -f "$output" ]] && [[ $(stat -c%s "$output" 2>/dev/null || stat -f%z "$output" 2>/dev/null) -gt 1000000 ]]; then
    local size=$(du -h "$output" 2>/dev/null | cut -f1)
    echo "  [跳过] hfi_${year}.tif 已存在 ($size)"
    return 0
  fi
  
  echo "  [下载] hfi_${year}.tif ..."
  echo "         URL: $url"
  
  if curl -L --max-time 1800 -o "$output" "$url" 2>/dev/null; then
    local size=$(du -h "$output" 2>/dev/null | cut -f1)
    echo "  [完成] hfi_${year}.tif ($size)"
    return 0
  else
    echo "  [失败] hfi_${year}.tif"
    rm -f "$output"
    return 1
  fi
}

# 2001-2014: HF3-1方法
echo ""
echo "[1/3] 下载 2001-2014 (HF3-1)..."
for year in $(seq 2001 2014); do
  url="https://storage.googleapis.com/hii-export/${year}-01-01/hii_${year}-01-01.tif"
  download_hfi "$year" "$url"
done

# 2015-2020: HF3-2方法 (推荐)
echo ""
echo "[2/3] 下载 2015-2020 (HF3-2, 含OSM数据)..."
for year in $(seq 2015 2020); do
  url="https://storage.googleapis.com/hii-export/${year}-01-01/hii_${year}-01-01.tif"
  download_hfi "$year" "$url"
done

# 2021-2024: 使用2020年数据作为近似（如果需要）
echo ""
echo "[3/3] 处理 2021-2024..."
echo "  注意: WCS HFI数据只到2020年。"
echo "  对于2021-2024年，建议:"
echo "    1. 使用2020年数据作为静态近似"
echo "    2. 或使用其他人类活动数据 (如GHSL Built-up)"

# 可选：创建符号链接
echo ""
read -p "是否创建2021-2024的符号链接指向2020年数据? (y/N): " create_links
if [[ "$create_links" =~ ^[Yy]$ ]]; then
  for year in $(seq 2021 2024); do
    if [[ ! -f "$DOWNLOAD_DIR/hfi_${year}.tif" ]]; then
      ln -sf "$DOWNLOAD_DIR/hfi_2020.tif" "$DOWNLOAD_DIR/hfi_${year}.tif"
      echo "  [链接] hfi_${year}.tif -> hfi_2020.tif"
    fi
  done
fi

echo ""
echo "========================================"
echo "下载完成！"
echo "========================================"
echo "文件列表:"
ls -lh "$DOWNLOAD_DIR"/hfi_*.tif 2>/dev/null || echo "  无tif文件"

TOTAL=$(ls "$DOWNLOAD_DIR"/hfi_*.tif 2>/dev/null | wc -l)
echo ""
echo "总计: $TOTAL 个tif文件"

# 统计缺失的年份
echo ""
echo "数据覆盖情况:"
for year in $(seq 2000 2024); do
  if [[ -f "$DOWNLOAD_DIR/hfi_${year}.tif" ]]; then
    size=$(du -h "$DOWNLOAD_DIR/hfi_${year}.tif" 2>/dev/null | cut -f1)
    echo "  ✓ $year ($size)"
  else
    echo "  ✗ $year (缺失)"
  fi
done

echo ""
echo "========================================"
echo "重要说明:"
echo "1. 2000年数据在WCS数据集中缺失"
echo "   - 建议: 使用2001年数据代替"
echo "   - 或从Figshare获取 (Mu et al. 2022)"
echo ""
echo "2. 2019-2024年数据说明:"
echo "   - 2019-2020: WCS提供"
echo "   - 2021-2024: 使用2020年静态数据"
echo "   - 如需动态数据，建议使用GHSL Built-up"
echo ""
echo "3. 数据格式: Cloud-Optimized GeoTIFF"
echo "   - 投影: WGS84 (EPSG:4326)"
echo "   - 分辨率: 约1km"
echo "   - 值域: 0-50+ (人类影响指数)"
echo "========================================"
