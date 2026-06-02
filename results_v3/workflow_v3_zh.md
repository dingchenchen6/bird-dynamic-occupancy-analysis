# v3 Workflow Document — v3_full_200sp_ar1_spatial_extended

## 项目概述
- **研究主题**：2000-2024年中国鸟类群落时空动态，stMsPGOcc 动态占有率模型校正检测偏差
- **核心模型**：spOccupancy::stMsPGOcc() + AR1 时间随机效应 + NNGP 空间随机效应（exponential covariance）
- **数据规模**：N/A 条去重记录，200 种建模物种，1247 个100km网格
- **分析期数**：N/A 个 5 年 primary period

## Pipeline 阶段

### Stage 1: 数据合并与去重（01_merge_birdwatch_ebird.R）
- 跨源去重（中国观鸟记录中心 + eBird）
- source-aware hash 去重策略
- 年份范围：2000-2024

### Stage 2: 调查历史构建（02_build_survey_history.R）
- 繁殖季过滤（4-5-6-7-8 月）
- 5 年 primary period 划分
- 检测协变量：log_events + has_duration

### Stage 3: 环境变量准备（03_prepare_environment.R）
- 4 组驱动变量：气候、地形+栖息地、人为、空间
- 标准化至 z-score

### Stage 3b: 性状扩展（03b_extend_traits.R）
- 基础 6 性状（body_mass, clutch_size, longevity, maturity, HWI, range_size）
- **新增** diet_specialization（Morelli et al. 2021，EltonTraits Shannon 逆指数）
- **新增** habitat_breadth（IUCN Red List 栖息地数量）
- 同物异名匹配 + 随机森林插值

### Stage 4: stMsPGOcc 模型拟合（04_run_stMsPGOcc_main.R）
- 空间多物种动态占有率模型（stMsPGOcc）
- NNGP 近似（N.neighbors=5，cov.model=exponential）
- 4 链 MCMC（batch=150, burn=5000, thin=2）
- 收敛：R-hat < 1.05: N/A%，ESS > 200: N/A%

### Stage 5: 后处理多样性（05_postprocess_diversity.R）
- 校正丰富度：80.35 → 102.74
- brms 驱动因子回归（gp 空间高斯过程）
- Baselga 时间 β 多样性分解

### Stage 6: 出版级图表（06/06b）
- Nature/Science 风格图表（Arial 7.5pt，89/120/183mm）
- varpart 纯分数：Climate=N/A%，Topo+Habitat=N/A%，Human=N/A%，Space=N/A%
- RF 排列重要性（100 draws × 1000 trees）：Top-3 = N/A

## 五个科学问题
1. **Q1** 多样性趋势：校正丰富度、Shannon、PD、功能多样性
2. **Q2** 定殖 vs 灭绝：物种水平动态
3. **Q3** 环境驱动因子：varpart（线性）+ RF（非线性）互补
4. **Q4** 性状解释力：含 diet_specialization + habitat_breadth
5. **Q5** 校正 vs 朴素：占有率校正 vs 原始计数
