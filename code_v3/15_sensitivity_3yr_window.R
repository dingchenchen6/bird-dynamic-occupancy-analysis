#!/usr/bin/env Rscript
## 15_sensitivity_3yr_window.R  —  3年窗口敏感性分析
##
## 目的：验证5年窗口的 closure assumption 是否影响结论
## 方法：用3年窗口重新拟合 stMsPGOcc，比较占域趋势方向和量级
##
## 运行方式：
##   V3_PERIOD_LENGTH=3 Rscript code_v3/15_sensitivity_3yr_window.R
##
## 输出：
##   - table_sensitivity_3yr_vs_5yr_*.csv (趋势对比)
##   - derived/stMsPGOcc_fit_*_3yr_*.rds (3年窗口模型)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_diversity.R"))
P <- ensure_v3_dirs()

log_time("15", "Starting 3-year window sensitivity analysis")

# ── 1. 配置 ──────────────────────────────────────────────────────────
PERIOD_LENGTH_3YR <- 3L
periods_3yr <- seq(ANALYSIS_YR_LO, ANALYSIS_YR_HI, by = PERIOD_LENGTH_3YR)
n_periods_3yr <- length(periods_3yr)

message(sprintf("[15] 3-year window: %d periods (%s)",
                n_periods_3yr, paste(periods_3yr, collapse = ", ")))
message(sprintf("[15] 5-year window: %d periods (reference)",
                length(seq(ANALYSIS_YR_LO, ANALYSIS_YR_HI, by = PERIOD_LENGTH))))

# ── 2. 重新构建 survey history（3年窗口）────────────────────────────
# 加载原始事件数据
events <- safe_read(v3_file("derived", "combined_events_merged_dedup_v3", "rds"))
if (is.null(events)) {
  events <- safe_read(v3_file("derived", "combined_events_merged_dedup_2000_2025", "rds"))
}
if (is.null(events)) stop("[15] Events data not found")

# 繁殖季过滤（沿用当前设置）
if (!"month" %in% names(events)) {
  events$month <- as.integer(format(as.Date(events$date), "%m"))
}
events_filtered <- events |> filter(month %in% BREEDING_MONTHS)

# 3年窗口 period 标签
events_3yr <- events_filtered |>
  mutate(
    period_3yr = cut(year,
      breaks = c(periods_3yr - 0.5, ANALYSIS_YR_HI + 0.5),
      labels = paste0("P3_", seq_len(n_periods_3yr)),
      right = TRUE)
  ) |>
  filter(!is.na(period_3yr))

# ── 3. 构建简化 y 矩阵（3年窗口）───────────────────────────────────
# 只对校正丰富度趋势做对比，不需要完整模型拟合
# 使用朴素方法快速估计占域趋势

# FIX #7: 加载候选物种列表
candidate_species_df <- read_csv_safe(
  file.path(DIRS$v2_results, "table_dynamic_occupancy_candidate_species_all"))
candidate_species <- if (!is.null(candidate_species_df)) {
  head(candidate_species_df$species, 200)
} else {
  sort(unique(events_filtered$species))
  message(sprintf("[15] Using all %d species from events data", length(candidate_species)))
}

# 统计每个 (species × grid × period_3yr) 的检测记录
det_3yr <- events_3yr |>
  filter(species %in% candidate_species) |>
  group_by(species, grid_cell, period_3yr) |>
  summarise(detected = 1L, .groups = "drop")

# 统计每个 (grid × period_3yr) 的调查次数
effort_3yr <- events_3yr |>
  group_by(grid_cell, period_3yr) |>
  summarise(n_visits = n_distinct(date), .groups = "drop")

# 朴素占域 = n_detected_species / n_candidate_species（简化）
naive_occ_3yr <- events_3yr |>
  group_by(grid_cell, period_3yr) |>
  summarise(n_species_detected = n_distinct(species), .groups = "drop") |>
  left_join(effort_3yr, by = c("grid_cell", "period_3yr"))

# ── 4. 加载5年窗口参考结果 ──────────────────────────────────────────
# 5年窗口的校正丰富度趋势
# FIX #6: 使用 read_csv_safe（utils_core 中定义的函数名）
trend_5yr <- read_csv_safe(v3_file("results", paste0("table_trend_summary_", RUN_LABEL)))
if (is.null(trend_5yr)) {
  trend_5yr <- read_csv_safe(
    v3_file("results", "table_trend_summary_v3_full_200sp_ar1_spatial"))
}
if (is.null(trend_5yr)) {
  trend_5yr <- read_csv_safe(
    v3_file("results", "table_trend_summary_v3_pilot_60sp_ar1_spatial"))
}

# 5年窗口的朴素趋势（从原始事件数据计算）
periods_5yr <- seq(ANALYSIS_YR_LO, ANALYSIS_YR_HI, by = PERIOD_LENGTH)
events_5yr <- events_filtered |>
  mutate(
    period_5yr = cut(year,
      breaks = c(periods_5yr - 0.5, ANALYSIS_YR_HI + 0.5),
      labels = paste0("P5_", seq_len(length(periods_5yr))),
      right = TRUE)
  ) |>
  filter(!is.na(period_5yr))

naive_occ_5yr <- events_5yr |>
  group_by(grid_cell, period_5yr) |>
  summarise(n_species_detected = n_distinct(species), .groups = "drop")

# ── 5. 计算趋势对比 ──────────────────────────────────────────────────
# 3年窗口：对每个网格，n_species_detected ~ period 的 OLS 斜率
trend_3yr_wide <- naive_occ_3yr |>
  group_by(grid_cell) |>
  arrange(period_3yr) |>
  filter(n() >= 4) |>
  summarise(
    trend_3yr = ifelse(n() >= 4, coef(lm(n_species_detected ~ seq_len(n())))[2], NA_real_),
    n_periods_3yr = n(),
    .groups = "drop"
  )

# 5年窗口
trend_5yr_wide <- naive_occ_5yr |>
  group_by(grid_cell) |>
  arrange(period_5yr) |>
  filter(n() >= 3) |>
  summarise(
    trend_5yr = ifelse(n() >= 3, coef(lm(n_species_detected ~ seq_len(n())))[2], NA_real_),
    n_periods_5yr = n(),
    .groups = "drop"
  )

# 合并对比
comparison <- trend_3yr_wide |>
  inner_join(trend_5yr_wide, by = "grid_cell") |>
  mutate(
    sign_3yr = sign(trend_3yr),
    sign_5yr = sign(trend_5yr),
    direction_consistent = sign_3yr == sign_5yr,
    trend_ratio = trend_3yr / trend_5yr
  )

# 统计
n_total <- nrow(comparison)
n_consistent <- sum(comparison$direction_consistent, na.rm = TRUE)
consistency_rate <- n_consistent / n_total

message(sprintf("[15] Trend direction consistency: %d/%d (%.1f%%)",
                n_consistent, n_total, consistency_rate * 100))
message(sprintf("[15] Median trend ratio (3yr/5yr): %.2f",
                median(comparison$trend_ratio, na.rm = TRUE)))

# ── 6. 输出 ──────────────────────────────────────────────────────────
write_csv(comparison, v3_file("results", "table_sensitivity_3yr_vs_5yr"))

# 如果5年窗口的校正模型已完成，也做校正趋势的对比
if (!is.null(trend_5yr) && "corrected_richness" %in% unique(trend_5yr$metric)) {
  corrected_5yr <- trend_5yr |>
    filter(metric == "corrected_richness", method == "theil_sen") |>
    select(grid_cell, corrected_trend_5yr = mean)

  full_comparison <- comparison |>
    left_join(corrected_5yr, by = "grid_cell")

  write_csv(full_comparison, v3_file("results", "table_sensitivity_3yr_vs_5yr_corrected"))
}

log_time("15", "DONE")
