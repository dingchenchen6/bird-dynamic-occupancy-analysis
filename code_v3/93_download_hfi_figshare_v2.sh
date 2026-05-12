#!/usr/bin/env bash
## 93_download_hfi_figshare_v2.sh  —  从Figshare下载HFI 2000-2018数据 (改进版)
##
## 数据集: Mu et al. 2022, Scientific Data
## DOI: 10.6084/m9.figshare.16571064.v5
##
## 注意: HFI数据只到2018年，2019-2024需要:
##   1. 使用2018年数据作为近似
##   2. 或使用其他人类活动数据源 (如GHSL, GHS-BUILT)

set -euo pipefail

DOWNLOAD_DIR="${1:-/home/dingchenchen/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3/data/hfi_figshare}"
mkdir -p "$DOWNLOAD_DIR"

echo "========================================"
echo "HFI Figshare 下载脚本 v2"
echo "目标目录: $DOWNLOAD_DIR"
echo "========================================"

# Figshare API endpoint
ARTICLE_ID="16571064"
API_URL="https://api.figshare.com/v2/articles/${ARTICLE_ID}/files"

echo "[1/2] 获取文件列表 from Figshare API..."
FILES_JSON=$(curl -s "$API_URL")

# 检查API响应
if [[ -z "$FILES_JSON" ]] || [[ "$FILES_JSON" == "null" ]]; then
  echo "ERROR: 无法获取Figshare文件列表"
  echo "尝试备用方法..."
  
  # 备用：直接构造URL（已知文件ID范围）
  # Figshare文件ID通常在 30000000-40000000 范围
  echo "使用备用下载方法..."
  
  # 已知的一些文件ID（需要验证）
  # 通过网页检查获取真实文件ID
  for year in $(seq 2000 2018); do
    output="$DOWNLOAD_DIR/hfp${year}.zip"
    [[ -f "$output" ]] && continue
    
    # 尝试多个可能的URL模式
    urls=(
      "https://figshare.com/ndownloader/files/30433885/hfp${year}.zip"
      "https://ndownloader.figshare.com/files/30433885/hfp${year}.zip"
    )
    
    for url in "${urls[@]}"; do
      echo "  尝试: $url"
      if curl -L --max-time 60 -o "$output" "$url" 2>/dev/null; then
        if [[ -s "$output" ]] && [[ $(stat -c%s "$output" 2>/dev/null || stat -f%z "$output" 2>/dev/null) -gt 1000 ]]; then
          size=$(du -h "$output" 2>/dev/null | cut -f1)
          echo "  [成功] hfp${year}.zip ($size)"
          break
        else
          rm -f "$output"
        fi
      fi
    done
  done
else
  echo "[2/2] 解析并下载文件..."
  
  # 使用Python解析JSON
  python3 << 'PYEOF'
import sys, json, os, subprocess

download_dir = os.environ.get('DOWNLOAD_DIR', '/tmp')
files_json = """${FILES_JSON}"""

try:
    files = json.loads(files_json)
    print(f"找到 {len(files)} 个文件")
    
    for f in files:
        name = f.get('name', '')
        url = f.get('download_url', '')
        
        if not name.endswith('.zip'):
            continue
            
        output = os.path.join(download_dir, name)
        if os.path.exists(output) and os.path.getsize(output) > 1000:
            print(f"  [跳过] {name} 已存在")
            continue
        
        print(f"  [下载] {name} ...")
        try:
            result = subprocess.run(
                ['curl', '-L', '--max-time', '600', '-o', output, url],
                capture_output=True, text=True, timeout=630
            )
            if result.returncode == 0 and os.path.exists(output):
                size = os.path.getsize(output)
                if size > 1000:
                    size_mb = size / (1024*1024)
                    print(f"  [完成] {name} ({size_mb:.1f} MB)")
                else:
                    print(f"  [失败] {name} 文件太小")
                    os.remove(output)
            else:
                print(f"  [失败] {name}")
                if os.path.exists(output):
                    os.remove(output)
        except Exception as e:
            print(f"  [错误] {name}: {e}")
            if os.path.exists(output):
                os.remove(output)
                
except Exception as e:
    print(f"ERROR parsing JSON: {e}")
PYEOF
fi

echo ""
echo "========================================"
echo "下载完成！"
echo "========================================"
ls -lh "$DOWNLOAD_DIR"/*.zip 2>/dev/null || echo "无zip文件"
TOTAL=$(ls "$DOWNLOAD_DIR"/*.zip 2>/dev/null | wc -l)
echo "总计: $TOTAL/19 个文件"

# 检查缺失的年份
echo ""
echo "缺失的年份:"
for year in $(seq 2000 2018); do
  if [[ ! -f "$DOWNLOAD_DIR/hfp${year}.zip" ]]; then
    echo "  - $year"
  fi
done

# 创建符号链接：2019-2024指向2018（如果需要）
echo ""
echo "注意: HFI数据只到2018年。"
echo "对于2019-2024年，建议:"
echo "  1. 使用2018年数据作为静态近似"
echo "  2. 或使用GHSL Built-up数据作为替代"
