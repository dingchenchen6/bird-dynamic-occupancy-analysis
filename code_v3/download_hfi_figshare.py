#!/usr/bin/env python3
"""
download_hfi_figshare.py - 从Figshare下载HFI 2000-2018数据

数据集: Mu et al. 2022, Scientific Data
DOI: 10.6084/m9.figshare.16571064
"""

import os
import sys
import json
import time
import urllib.request
import urllib.error

# 配置
DOWNLOAD_DIR = os.environ.get("HFI_DOWNLOAD_DIR", "/home/dingchenchen/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3/data/hfi_figshare")
ARTICLE_ID = "16571064"

os.makedirs(DOWNLOAD_DIR, exist_ok=True)

# 尝试通过Figshare API获取文件列表
def get_files_from_api():
    """通过Figshare API获取文件列表"""
    api_url = f"https://api.figshare.com/v2/articles/{ARTICLE_ID}/files"
    
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
    }
    
    try:
        req = urllib.request.Request(api_url, headers=headers)
        with urllib.request.urlopen(req, timeout=30) as response:
            data = json.loads(response.read().decode('utf-8'))
            return data
    except Exception as e:
        print(f"API请求失败: {e}")
        return None

def download_file(url, output_path, max_retries=3):
    """下载单个文件"""
    if os.path.exists(output_path) and os.path.getsize(output_path) > 1000:
        size_mb = os.path.getsize(output_path) / (1024 * 1024)
        print(f"  [跳过] {os.path.basename(output_path)} 已存在 ({size_mb:.1f} MB)")
        return True
    
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
    }
    
    for attempt in range(max_retries):
        try:
            print(f"  [下载] {os.path.basename(output_path)} ... (尝试 {attempt + 1}/{max_retries})")
            req = urllib.request.Request(url, headers=headers)
            
            with urllib.request.urlopen(req, timeout=600) as response:
                with open(output_path, 'wb') as f:
                    while True:
                        chunk = response.read(1024 * 1024)  # 1MB chunks
                        if not chunk:
                            break
                        f.write(chunk)
            
            size_mb = os.path.getsize(output_path) / (1024 * 1024)
            print(f"  [完成] {os.path.basename(output_path)} ({size_mb:.1f} MB)")
            return True
            
        except Exception as e:
            print(f"  [失败] 尝试 {attempt + 1}: {e}")
            if os.path.exists(output_path):
                os.remove(output_path)
            time.sleep(5)
    
    return False

def main():
    print("=" * 60)
    print("HFI Figshare 下载脚本")
    print(f"目标目录: {DOWNLOAD_DIR}")
    print("=" * 60)
    
    # 方法1: 通过API获取文件列表
    print("\n[1/2] 尝试通过Figshare API获取文件列表...")
    files = get_files_from_api()
    
    if files:
        print(f"找到 {len(files)} 个文件")
        
        success_count = 0
        for file_info in files:
            name = file_info.get('name', '')
            url = file_info.get('download_url', '')
            
            if not name.endswith('.zip'):
                continue
            
            output_path = os.path.join(DOWNLOAD_DIR, name)
            if download_file(url, output_path):
                success_count += 1
            
            time.sleep(1)  # 礼貌延迟
        
        print(f"\n成功下载: {success_count} 个文件")
    else:
        print("API获取失败，尝试备用方法...")
        
        # 方法2: 尝试已知的直接URL模式
        # Figshare文件ID需要通过其他方式获取
        print("\n[2/2] 尝试通过备用URL下载...")
        
        # 已知的文件ID（需要从网页获取）
        # 这里使用一个示例，实际ID需要通过浏览器检查获取
        known_files = {
            # 年份: 文件ID (需要更新)
        }
        
        if not known_files:
            print("错误: 没有已知的文件ID。请手动从Figshare网页获取下载链接。")
            print("\n手动下载步骤:")
            print("1. 访问 https://figshare.com/articles/figure/An_annual_global_terrestrial_Human_Footprint_dataset_from_2000_to_2018/16571064")
            print("2. 点击每个年份的'Download'按钮")
            print("3. 或使用浏览器开发者工具获取直接下载链接")
            sys.exit(1)

if __name__ == "__main__":
    main()
