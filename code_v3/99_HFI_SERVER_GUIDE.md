# HFI数据下载与服务器运行指南

## 概述

本指南说明如何在服务器(162.105.149.23)上下载HFI数据并运行v3 pipeline。

## 数据源选择

### 方案1: WCS HFI (推荐，可直接服务器下载)
- **来源**: Wildlife Conservation Society
- **网址**: https://www.wcshumanfootprint.org/data-access
- **时间范围**: 2001-2020
- **格式**: Cloud-Optimized GeoTIFF
- **投影**: WGS84 (EPSG:4326)
- **分辨率**: ~1km
- **优点**: 可直接在服务器下载，无需手动操作
- **缺点**: 2000年缺失，2021-2024需用2020年代替

### 方案2: Figshare/Mu et al. (原始论文数据)
- **来源**: Mu et al. 2022, Scientific Data
- **DOI**: 10.6084/m9.figshare.16571064
- **时间范围**: 2000-2018
- **格式**: zip (内含GeoTIFF)
- **投影**: Mollweide等积投影
- **分辨率**: 1km
- **优点**: 有2000年数据，与论文一致
- **缺点**: Figshare反爬虫，需手动获取下载链接

### 方案3: 百度网盘 (中国农大UEMM团队)
- **链接**: https://pan.baidu.com/s/1X9OB1xjQdo2x4uk6lHJHNg
- **密码**: tvm7
- **时间范围**: 2000-2020
- **优点**: 时间范围完整
- **缺点**: 服务器无法直接访问，需本机下载后上传

## 快速开始

### 方法A: 全自动 (推荐)

在服务器上执行:

```bash
# 1. 设置环境变量
export V3_GRID_SIZE_KM="10"
export V3_MIN_VISITS="1"
export V3_MIN_SPECIES_GRID="3"
export V3_N_OMP_THREADS="16"

# 2. 创建tmux会话并启动
cd ~/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3
bash code_v3/96_tmux_launch_server.sh

# 3. 分离会话 (按 Ctrl+B, 然后按 D)
# 现在可以安全断开SSH，任务会在后台继续
```

### 方法B: 分步执行

```bash
# 1. 下载HFI数据
bash code_v3/93_download_hfi_wcs.sh \
  ~/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3/data/hfi_wcs

# 2. 运行完整pipeline
export V3_HFI_WCS_DIR=~/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3/data/hfi_wcs
bash code_v3/95_server_hfi_and_pipeline.sh
```

## 监控进度

### 查看tmux会话
```bash
tmux ls                                    # 列出所有会话
tmux attach -t hfi_download               # 查看HFI下载进度
tmux attach -t bird_pipeline              # 查看Pipeline进度
tmux detach                               # 分离会话 (Ctrl+B, D)
```

### 查看日志
```bash
# HFI下载日志
tail -f ~/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3/logs_v3/hfi_download_*.log

# Pipeline日志
tail -f ~/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3/logs_v3/pipeline_*.log

# 模型运行日志
tail -f ~/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3/logs_v3/04_run_stMsPGOcc_main_*.log
```

## 文件说明

| 脚本 | 用途 |
|------|------|
| `93_download_hfi_wcs.sh` | 从WCS下载HFI 2001-2020 |
| `93_download_hfi_figshare.sh` | Figshare下载脚本(自动) |
| `94_download_hfi_figshare_manual.sh` | Figshare下载辅助(手动) |
| `95_server_hfi_and_pipeline.sh` | HFI下载+Pipeline整合 |
| `96_tmux_launch_server.sh` | Tmux后台启动 |

## 注意事项

1. **2000年数据**: WCS数据集从2001年开始，如需2000年数据，建议:
   - 使用Figshare/Mu et al.数据
   - 或用2001年代替

2. **2019-2024年数据**: WCS提供到2020年，2021-2024建议:
   - 使用2020年静态数据
   - 或使用GHSL Built-up数据

3. **存储空间**: 每个tif文件约3-3.5GB，20年约需70GB

4. **网络稳定**: 使用tmux确保下载不受SSH断开影响

## 自动化进度检查

已设置自动化任务，每2小时检查一次进度:
- 任务ID: automation-1777921816720
- 检查内容: HFI下载进度、Pipeline阶段、模型进度、系统资源
