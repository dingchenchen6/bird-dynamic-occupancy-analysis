#!/usr/bin/env bash
## 94_download_hfi_figshare_manual.sh  —  Figshare HFI数据下载辅助脚本
##
## 由于Figshare API限制，此脚本提供手动下载指导
##
## 数据集: Mu et al. 2022, Scientific Data
## DOI: 10.6084/m9.figshare.16571064
## 时间范围: 2000-2018 (19个zip文件)

set -euo pipefail

DOWNLOAD_DIR="${1:-/home/dingchenchen/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3/data/hfi_figshare}"
mkdir -p "$DOWNLOAD_DIR"

echo "========================================"
echo "Figshare HFI 数据下载辅助脚本"
echo "========================================"
echo ""
echo "由于Figshare的反爬虫限制，自动下载可能失败。"
echo "请按以下步骤手动获取下载链接:"
echo ""
echo "步骤1: 访问Figshare文章页面"
echo "  URL: https://figshare.com/articles/figure/An_annual_global_terrestrial_Human_Footprint_dataset_from_2000_to_2018/16571064"
echo ""
echo "步骤2: 获取文件下载链接"
echo "  方法A - 浏览器开发者工具:"
echo "    1. 在浏览器中打开上述URL"
echo "    2. 按F12打开开发者工具 -> Network标签"
echo "    3. 点击任意年份的'Download'按钮"
echo "    4. 在Network中找到以'ndownloader'开头的请求"
echo "    5. 右键 -> Copy -> Copy link address"
echo ""
echo "  方法B - 使用wget/curl递归下载:"
echo "    wget --recursive --no-parent --accept '*.zip' \\"
echo "         https://figshare.com/articles/figure/An_annual_global_terrestrial_Human_Footprint_dataset_from_2000_to_2018/16571064"
echo ""
echo "步骤3: 将下载的zip文件放入目录:"
echo "  $DOWNLOAD_DIR"
echo ""
echo "========================================"
echo "备用数据源 (推荐直接在服务器下载):"
echo "========================================"
echo ""
echo "1. WCS HFI (2001-2020)"
echo "   运行: bash code_v3/93_download_hfi_wcs.sh"
echo "   优点: 可直接下载，无需手动操作"
echo "   缺点: 2000年缺失，2019-2020为近似"
echo ""
echo "2. 百度网盘 (2000-2018/2020)"
echo "   链接: https://pan.baidu.com/s/1X9OB1xjQdo2x4uk6lHJHNg"
echo "   密码: tvm7"
echo "   注意: 服务器无法直接访问百度网盘，需先下载到本机再上传"
echo ""
echo "3. 中国农业大学UEMM团队"
echo "   页面: https://www.x-mol.com/groups/li_xuecao/news/48145"
echo ""
echo "========================================"
echo "数据对比:"
echo "========================================"
echo ""
echo "| 特性 | Figshare (Mu) | WCS | 百度网盘 |"
echo "|------|---------------|-----|----------|"
echo "| 时间范围 | 2000-2018 | 2001-2020 | 2000-2020 |"
echo "| 格式 | zip (tif) | tif | 未知 |"
echo "| 投影 | Mollweide | WGS84 | Mollweide |"
echo "| 分辨率 | 1km | 1km | 1km |"
echo "| 服务器下载 | 困难 | 容易 | 不可 |"
echo ""
echo "========================================"

# 检查已有文件
echo ""
echo "当前目录已有文件:"
ls -lh "$DOWNLOAD_DIR"/*.zip 2>/dev/null || echo "  无zip文件"

# 如果用户提供了下载链接文件，尝试下载
LINKS_FILE="$DOWNLOAD_DIR/download_links.txt"
if [[ -f "$LINKS_FILE" ]]; then
  echo ""
  echo "发现下载链接文件: $LINKS_FILE"
  echo "开始批量下载..."
  
  while IFS= read -r line; do
    [[ -z "$line" ]] && continue
    [[ "$line" =~ ^# ]] && continue
    
    # 格式: 年份 URL
    year=$(echo "$line" | awk '{print $1}')
    url=$(echo "$line" | awk '{print $2}')
    
    output="$DOWNLOAD_DIR/hfp${year}.zip"
    if [[ -f "$output" ]]; then
      echo "  [跳过] hfp${year}.zip 已存在"
      continue
    fi
    
    echo "  [下载] hfp${year}.zip ..."
    if curl -L --max-time 600 -o "$output" "$url" 2>/dev/null; then
      size=$(du -h "$output" 2>/dev/null | cut -f1)
      echo "  [完成] hfp${year}.zip ($size)"
    else
      echo "  [失败] hfp${year}.zip"
      rm -f "$output"
    fi
    
    sleep 2
  done < "$LINKS_FILE"
  
  echo ""
  echo "批量下载完成！"
fi

echo ""
echo "========================================"
echo "提示: 如需帮助获取Figshare下载链接，"
echo "      请提供服务器SSH访问方式，"
echo "      我可以远程协助操作浏览器获取链接。"
echo "========================================"
