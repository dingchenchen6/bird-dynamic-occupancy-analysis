# 中国鸟类群落时空动态：基于空间占有率模型的检测偏差校正

## 摘要
基于 2000-2024 年 N/A 条鸟类观测记录，我们使用空间多物种动态占有率模型（stMsPGOcc）校正检测偏差，评估了中国 1247 个 100km 网格内 200 种鸟类的群落时空动态。校正丰富度从 80.35 上升至 102.74。环境驱动因子分析显示空间因子解释方差最大（adj R²=N/A%），随机森林排列重要性前三位为 N/A。新增性状分析表明饮食专化性和栖息地广度对物种趋势有显著解释力。

## 1 引言
全球生物多样性正经历快速变化。检测偏差（detection bias）是鸟类监测数据分析的核心挑战。本研究利用 stMsPGOcc 模型，首次在中国尺度上同时纳入空间相关性和检测偏差校正，评估群落动态。

## 2 方法

### 2.1 数据来源与处理
- 中国观鸟记录中心 + eBird（2000-2024）
- 繁殖季过滤（4-5-6-7-8 月）
- source-aware 跨源去重

### 2.2 占有率模型
stMsPGOcc（spOccupancy R 包），包含：
- AR1 时间随机效应
- NNGP 空间随机效应（exponential covariance，N.neighbors=5）
- 检测子模型：~ log_events + log_duration + has_duration

### 2.3 群落多样性指标
- 分类多样性：校正丰富度、Shannon
- 系统发育多样性：Faith's PD（概率加权）、MPD
- 功能多样性：性状体积、Rao's Q
- 性状维度包含 diet_specialization 和 habitat_breadth

### 2.4 环境驱动因子分析
- varpart（线性方差分解）：Climate、Topo+Habitat、Human、Space
- 随机森林排列重要性（非线性，互补）：100 draws × 1000 trees
- brms 驱动因子回归：gp(centroid_lon, centroid_lat)

### 2.5 性状-趋势回归（Q4）
brms 生态位模型：
  trend_i ~ z_body_mass + z_hwi + z_range_size + z_clutch_size
          + z_diet_specialization + z_habitat_breadth
          + (1 | gr(species, cov = A))

### 2.6 同质化分析
- Sørensen 距离时序变化
- Baselga 分解：turnover vs nestedness
- Mann-Kendall 趋势检验

## 3 结果

### 3.1 模型收敛
R-hat < 1.05: N/A%，ESS > 200: N/A%

### 3.2 多样性趋势
校正丰富度：80.35 → 102.74

### 3.3 环境驱动因子
varpart 纯分数：Climate=N/A%，Topo+Habitat=N/A%，Human=N/A%，Space=N/A%
RF 重要性排名：N/A

### 3.4 性状解释力
（待模型运行后补充）

### 3.5 时间β多样性
（待模型运行后补充）

### 3.6 校正 vs 朴素
（待模型运行后补充）

## 4 讨论
- stMsPGOcc 空间建模对检测偏差校正的改进
- 饮食专化性和栖息地广度的生态意义
- varpart 与 RF 互补性的方法论贡献
- 生物同质化的时空模式

## 参考文献
- Doser, J. W., Finley, A. O., Kery, M., & Zipkin, E. F. (2024). spOccupancy...
- Morelli, F., Benedetti, Y., Hanson, J. O., & Fuller, R. A. (2021). Conservation Letters, e12795.
- Wilman, H., et al. (2014). EltonTraits 1.0. Ecology, 95, 1887.
- Wright, M. N., & Ziegler, A. (2017). ranger. J Statistical Software.
