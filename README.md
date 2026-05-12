# 中国鸟类群落多维时空动态研究

## Multidimensional Spatiotemporal Dynamics of Chinese Avian Communities

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://cran.r-project.org/)
[![spOccupancy](https://img.shields.io/badge/spOccupancy-v0.8.0-green.svg)](https://doserlab.com/packages/spoccupancy/)
[![brms](https://img.shields.io/badge/brms-v2.20+-green.svg)](https://paul-buerkner.github.io/brms/)

---

## 项目概述 / Project Overview

本项目构建首个整合**时空多物种动态占用模型（stMsPGOcc）**与**多维多样性后验概率加权估计**的贝叶斯分析框架，系统量化2000–2024年间中国鸟类群落在分类学、功能和系统发育三个维度上的时空动态格局，识别功能同质化风险热点，归因气候、土地利用和人类干扰的独立与交互效应。

This project establishes the first Bayesian analytical framework integrating **spatiotemporal multi-species dynamic occupancy models (stMsPGOcc)** with **posterior probability-weighted estimation of multidimensional diversity**, to systematically quantify the spatiotemporal dynamics of Chinese avian communities across taxonomic, functional, and phylogenetic dimensions from 2000–2024.

---

## 核心创新 / Key Innovations

1. **检测概率校正**：使用 stMsPGOcc（AR1 + spatial NNGP）校正不完全检测，解决传统方法对稀有物种丰富度的系统性低估
2. **多维后验传播**：将功能多样性（FEve, FDiv, Rao's Q）和系统发育多样性（Faith's PD, MPD）与占用模型后验耦合
3. **概率加权β分解**：Baselga turnover/nestedness 分解纳入后验不确定性
4. **物种水平趋势系统分析**：200物种的扩张/稳定/收缩分类与性状-环境关联

---

## 分析流程 / Analysis Pipeline

```
Stage 01  数据整合      Birdwatch × eBird 融合
Stage 02  调查历史      多时期检测-非检测矩阵构建
Stage 03  环境数据      气候/地形/土地利用/人类干扰准备
Stage 04  核心模型      stMsPGOcc 拟合 (200 sp × 1247 sites × 5 periods)
Stage 05  后处理        多维多样性计算 + 趋势估计 + brms驱动回归
Stage 06  出版图集      Nature/Science 级别图集
Stage 07  稿件渲染      Markdown 稿件生成
Stage 10  PPTX          演示文稿
Stage 11  期刊DOCX      期刊格式文档
Stage 12  同质化地图    时空同质化风险地图
Stage 13  全球提案      政策提案文档
```

---

## 目录结构 / Repository Structure

```
bird_dynamic_occupancy_analysis/
├── code_v3/                          # v3 分析代码 (主版本)
│   ├── 00_config.R                   # 集中配置 (MCMC参数/路径/常量)
│   ├── 01_merge_birdwatch_ebird.R    # 数据整合
│   ├── 02_build_survey_history.R     # 调查历史构建
│   ├── 03_prepare_environment.R      # 环境数据准备
│   ├── 03b_extend_traits.R           # 性状数据扩展
│   ├── 03c_prepare_climate_change.R  # 气候变化数据
│   ├── 03d_prepare_landuse_change.R  # 土地利用变化
│   ├── 03e_prepare_hfi_change.R      # 人类足迹变化
│   ├── 04_run_stMsPGOcc_main.R       # 核心模型脚本
│   ├── 04a_run_stMsPGOcc_parallel.sh # 4链并行启动
│   ├── 04b_recover_diagnostics.R     # 诊断与合并
│   ├── 05_postprocess_diversity.R    # 后处理与多样性计算
│   ├── 06_figures_publication.R      # 出版级图集
│   ├── 06b_regenerate_driver_plots.R # 驱动因子图
│   ├── 07_render_manuscript.R        # 稿件渲染
│   ├── 08_post_full_run_pipeline.sh  # 完整后处理管线
│   ├── 09_extended_analyses.R        # 扩展分析
│   ├── 10_finalize_with_dashline_pptx.R  # PPTX生成
│   ├── 11_render_journal_docx.R      # 期刊DOCX
│   ├── 12_homogenization_spatiotemporal_maps.R  # 同质化地图
│   ├── 13_render_global_proposal_docx.R  # 全球提案
│   ├── 14_species_trait_regression.R     # 物种-性状回归
│   ├── 15_sensitivity_3yr_window.R       # 3年窗口敏感性
│   ├── 15b_sensitivity_breeding_season.R # 繁殖季敏感性
│   ├── utils_*.R                       # 工具函数库
│   └── 90_server_run_v3_pipeline.sh    # 服务器一键运行
├── .gitignore
└── README.md
```

---

## 依赖包 / Dependencies

### R 核心包
```r
# 统计建模
spOccupancy     # stMsPGOcc 时空占用模型
brms            # 贝叶斯空间回归
rstan / cmdstanr # MCMC后端

# 生态学分析
ape             # 系统发育分析
vegan           # 生态多样性

# 数据处理
dplyr, tidyr, readr, tibble, stringr
sf, terra, exactextractr  # 空间分析

# 可视化
ggplot2, patchwork, cowplot

# 文档
officer, flextable, rmarkdown
```

### 系统要求
- R ≥ 4.0
- Stan/CmdStan ≥ 2.26
- 内存：≥ 64GB（模型拟合建议 ≥ 256GB）
- 存储：≥ 500GB（原始数据 + 后验样本）

---

## 快速开始 / Quick Start

### 1. 克隆仓库
```bash
git clone https://github.com/dingchenchen6/bird-dynamic-occupancy-analysis.git
cd bird-dynamic-occupancy-analysis
```

### 2. 配置环境变量
```bash
export BIRD_PROJECT_ROOT=$(pwd)
export V3_CODE_DIR=$(pwd)/code_v3
```

### 3. 运行完整管线
```bash
# 本地运行
bash code_v3/98_smart_pipeline_05to13.sh

# 服务器运行（256核, 1TB内存）
bash code_v3/90_server_run_v3_pipeline.sh
```

### 4. 单独运行Stage
```bash
# 仅运行核心模型 (Stage 04)
Rscript code_v3/04_run_stMsPGOcc_main.R

# 仅运行后处理 (Stage 05)
Rscript code_v3/05_postprocess_diversity.R

# 仅运行图集 (Stage 06)
Rscript code_v3/06_figures_publication.R
```

---

## 模型配置 / Model Configuration

关键参数在 `code_v3/00_config.R` 中集中管理：

| 参数 | 值 | 说明 |
|------|-----|------|
| `FULL_N_BATCH` | 150 | MCMC batch数 |
| `FULL_N_BURN` | 5000 | burn-in迭代数 |
| `FULL_N_THIN` | 2 | thinning间隔 |
| `FULL_N_CHAINS` | 4 | 并行链数 |
| `PSI_MAX_DRAWS` | 400 | 后验抽取数 |
| `N_NEIGHBORS` | 5 | NNGP邻居数 |
| `COV_MODEL` | "exponential" | 空间协方差函数 |

---

## 数据说明 / Data

本项目使用的数据包括：

| 数据类型 | 来源 | 时间范围 |
|---------|------|---------|
| 鸟类检测记录 | Birdwatch + eBird | 2000–2024 |
| 气候数据 | WorldClim 2.1 | 2000–2024 |
| 土地利用 | MODIS MCD12Q1 | 2000–2024 |
| 人类足迹 | Venter et al. (2016) 更新版 | 2000–2024 |
| 地形 | SRTM/ASTER | 静态 |
| 性状 | EltonTraits 1.0 + HBW | 静态 |
| 系统发育 | Jetz et al. (2012) Hackett骨干 | 静态 |

**注意**：由于数据文件较大（>50GB），原始数据未包含在本仓库中。请联系作者获取数据访问权限。

---

## 作者 / Author

**丁晨晨 (Chenchen Ding)**
- 北京大学 城市与环境学院 生态研究中心
- Email: dingchenchen@pku.edu.cn
- GitHub: [@dingchenchen6](https://github.com/dingchenchen6)
- Website: [https://dingchenchen6.github.io/chenchen-ding/](https://dingchenchen6.github.io/chenchen-ding/)

---

## 许可 / License

本项目的分析代码遵循 [MIT License](LICENSE)。

数据使用遵循原始数据提供者的使用条款（eBird Basic Dataset Terms of Use, Birdwatch数据使用协议）。

---

## 引用 / Citation

如果您使用了本项目的代码或方法，请引用：

```bibtex
@software{ding2026birdoccupancy,
  author = {Ding, Chenchen},
  title = {Multidimensional Spatiotemporal Dynamics of Chinese Avian Communities},
  year = {2026},
  url = {https://github.com/dingchenchen6/bird-dynamic-occupancy-analysis}
}
```

---

## 致谢 / Acknowledgments

本项目得到以下支持和资助：
- 北京大学城市与环境学院
- 中国观鸟记录中心 (birdwatch.cn) 和 eBird 志愿者提供的长期监测数据
- spOccupancy 开发团队 (Doser et al., 2022)
