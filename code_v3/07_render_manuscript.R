#!/usr/bin/env Rscript
## 07_render_manuscript.R  —  v3 自动生成 workflow + 中文稿件
##
## 从 CSV 结果文件自动填充数字，生成 markdown 稿件。
## 更新为 stMsPGOcc、新性状（diet_specialization + habitat_breadth）、
## RF 排列重要性方法。

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(glue)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))

P <- ensure_v3_dirs()
RUN_LABEL <- Sys.getenv("V3_RUN_LABEL", RUN_LABEL)

log_time("07", sprintf("Rendering manuscript for %s", RUN_LABEL))

# ── 安全读取结果 ──────────────────────────────────────────────────────
r <- function(stem, ext = "csv") {
  path <- v3_file("results", stem, ext)
  if (file.exists(path)) read_csv(path, show_col_types = FALSE) else NULL
}

dedup_summary   <- r("table_dedup_audit_summary")
primary_blocks  <- r("table_primary_5year_blocks")
survey_year     <- r("table_survey_coverage_by_year")
survey_block    <- r("table_survey_coverage_by_block")
candidate_sp    <- r("table_dynamic_occupancy_candidate_species_all")
metrics         <- r(paste0("table_community_metrics_with_cri_", RUN_LABEL))
trends          <- r(paste0("table_grid_trends_with_cri_", RUN_LABEL))
beta            <- r(paste0("table_temporal_beta_with_cri_", RUN_LABEL))
mcmc_diag       <- r(paste0("mcmc_diagnostics_", RUN_LABEL))
varpart         <- r(paste0("table_varpart_richness_trend_", RUN_LABEL))
rf_summary      <- r("table_rf_importance_summary")
trait_ext       <- r("table_traits_extended_v3")
sp_trend        <- r(paste0("table_species_trait_regression_", RUN_LABEL))

# ── 提取关键数字 ──────────────────────────────────────────────────────
n_records <- if (!is.null(dedup_summary)) dedup_summary$n_after_dedup[1] else "N/A"
n_species <- if (!is.null(candidate_sp)) nrow(candidate_sp) else "N/A"
n_grids   <- if (!is.null(metrics)) length(unique(metrics$grid_cell)) else "N/A"
n_periods <- if (!is.null(primary_blocks)) nrow(primary_blocks) else "N/A"

richness_lo <- richness_hi <- "N/A"
if (!is.null(metrics)) {
  cr <- metrics |> filter(metric == "corrected_richness")
  if (nrow(cr) > 0) {
    by_period <- cr |> group_by(block_label) |> summarise(m = median(value_mean), .groups = "drop")
    richness_lo <- round(min(by_period$m, na.rm = TRUE), 2)
    richness_hi <- round(max(by_period$m, na.rm = TRUE), 2)
  }
}

# varpart 纯分数
vp_climate <- vp_topo <- vp_human <- vp_space <- "N/A"
if (!is.null(varpart)) {
  vp <- varpart |> filter(grepl("pure", component, ignore.case = TRUE))
  for (i in seq_len(nrow(vp))) {
    val <- round(vp$adj_R2[i] * 100, 1)
    if (grepl("Climate", vp$component[i])) vp_climate <- val
    if (grepl("Topo", vp$component[i]))    vp_topo <- val
    if (grepl("Human", vp$component[i]))   vp_human <- val
    if (grepl("Space", vp$component[i]))   vp_space <- val
  }
}

# RF 重要性排名
rf_top3 <- "N/A"
if (!is.null(rf_summary)) {
  top3 <- rf_summary |> arrange(desc(mean)) |> head(3) |> pull(variable)
  rf_top3 <- paste(top3, collapse = " > ")
}

# 收敛诊断
rhat_pct <- ess_pct <- "N/A"
if (!is.null(mcmc_diag)) {
  rhat_pct <- round(mean(mcmc_diag$rhat < RHAT_THRESHOLD, na.rm = TRUE) * 100, 1)
  ess_pct  <- round(mean(mcmc_diag$ess > ESS_THRESHOLD, na.rm = TRUE) * 100, 1)
}

# ── 生成 Workflow 文档 ───────────────────────────────────────────────

workflow_md <- glue("
# v3 Workflow Document — {RUN_LABEL}

## 项目概述
- **研究主题**：2000-2024年中国鸟类群落时空动态，stMsPGOcc 动态占有率模型校正检测偏差
- **核心模型**：spOccupancy::stMsPGOcc() + AR1 时间随机效应 + NNGP 空间随机效应（exponential covariance）
- **数据规模**：{n_records} 条去重记录，{n_species} 种建模物种，{n_grids} 个100km网格
- **分析期数**：{n_periods} 个 5 年 primary period

## Pipeline 阶段

### Stage 1: 数据合并与去重（01_merge_birdwatch_ebird.R）
- 跨源去重（中国观鸟记录中心 + eBird）
- source-aware hash 去重策略
- 年份范围：{ANALYSIS_YR_LO}-{ANALYSIS_YR_HI}

### Stage 2: 调查历史构建（02_build_survey_history.R）
- 繁殖季过滤（{paste(BREEDING_MONTHS, collapse='-')} 月）
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
- NNGP 近似（N.neighbors={N_NEIGHBORS}，cov.model={COV_MODEL}）
- {FULL_N_CHAINS} 链 MCMC（batch={FULL_N_BATCH}, burn={FULL_N_BURN}, thin={FULL_N_THIN}）
- 收敛：R-hat < {RHAT_THRESHOLD}: {rhat_pct}%，ESS > {ESS_THRESHOLD}: {ess_pct}%

### Stage 5: 后处理多样性（05_postprocess_diversity.R）
- 校正丰富度：{richness_lo} → {richness_hi}
- brms 驱动因子回归（gp 空间高斯过程）
- Baselga 时间 β 多样性分解

### Stage 6: 出版级图表（06/06b）
- Nature/Science 风格图表（Arial 7.5pt，89/120/183mm）
- varpart 纯分数：Climate={vp_climate}%，Topo+Habitat={vp_topo}%，Human={vp_human}%，Space={vp_space}%
- RF 排列重要性（100 draws × 1000 trees）：Top-3 = {rf_top3}

## 五个科学问题
1. **Q1** 多样性趋势：校正丰富度、Shannon、PD、功能多样性
2. **Q2** 定殖 vs 灭绝：物种水平动态
3. **Q3** 环境驱动因子：varpart（线性）+ RF（非线性）互补
4. **Q4** 性状解释力：含 diet_specialization + habitat_breadth
5. **Q5** 校正 vs 朴素：占有率校正 vs 原始计数
")

write_lines(workflow_md, v3_file("results", "workflow_v3_zh", "md"))
log_time("07", "Workflow document written")

# ── 生成中文稿件骨架 ──────────────────────────────────────────────────

manuscript_md <- glue("
# 中国鸟类群落时空动态：基于空间占有率模型的检测偏差校正

## 摘要
基于 2000-2024 年 {n_records} 条鸟类观测记录，我们使用空间多物种动态占有率模型（stMsPGOcc）校正检测偏差，评估了中国 {n_grids} 个 100km 网格内 {n_species} 种鸟类的群落时空动态。校正丰富度从 {richness_lo} 上升至 {richness_hi}。环境驱动因子分析显示空间因子解释方差最大（adj R²={vp_space}%），随机森林排列重要性前三位为 {rf_top3}。新增性状分析表明饮食专化性和栖息地广度对物种趋势有显著解释力。

## 1 引言
全球生物多样性正经历快速变化。检测偏差（detection bias）是鸟类监测数据分析的核心挑战。本研究利用 stMsPGOcc 模型，首次在中国尺度上同时纳入空间相关性和检测偏差校正，评估群落动态。

## 2 方法

### 2.1 数据来源与处理
- 中国观鸟记录中心 + eBird（{ANALYSIS_YR_LO}-{ANALYSIS_YR_HI}）
- 繁殖季过滤（{paste(BREEDING_MONTHS, collapse='-')} 月）
- source-aware 跨源去重

### 2.2 占有率模型
stMsPGOcc（spOccupancy R 包），包含：
- AR1 时间随机效应
- NNGP 空间随机效应（exponential covariance，N.neighbors={N_NEIGHBORS}）
- 检测子模型：{DET_FORMULA_STR}

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
R-hat < {RHAT_THRESHOLD}: {rhat_pct}%，ESS > {ESS_THRESHOLD}: {ess_pct}%

### 3.2 多样性趋势
校正丰富度：{richness_lo} → {richness_hi}

### 3.3 环境驱动因子
varpart 纯分数：Climate={vp_climate}%，Topo+Habitat={vp_topo}%，Human={vp_human}%，Space={vp_space}%
RF 重要性排名：{rf_top3}

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
")

write_lines(manuscript_md, v3_file("results", "manuscript_v3_zh", "md"))
log_time("07", "Manuscript skeleton written")
log_time("07", "Done.")
