#!/usr/bin/env Rscript
## 15b_sensitivity_breeding_season.R  —  繁殖季过滤敏感性分析
##
## 目的：评估繁殖季过滤（4-8月）vs 全年数据对占域估计的影响
##
## 核心问题：
##   1. 繁殖季过滤是否引入时间偏差（早期数据损失更多）？
##   2. 非繁殖季数据是否包含越冬/迁徙信号，混淆占域估计？
##   3. 过滤后 vs 全年的物种组成差异有多大？
##
## 输出：
##   - table_breeding_sensitivity_*.csv
##   - breeding_season_analysis_*.csv

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
P <- ensure_v3_dirs()

log_time("15b", "Starting breeding season sensitivity analysis")

# ── 1. 加载事件数据 ──────────────────────────────────────────────────
events_path <- v3_file("derived", "combined_events_dedup_v3", "rds")
if (!file.exists(events_path)) {
  events_path <- file.path(DIRS$v2_derived, "combined_events_merged_dedup_2000_2025.rds")
}
events <- safe_read(events_path)
if (is.null(events)) stop("[15b] Events data not found")

if (!"month" %in% names(events)) {
  events$month <- as.integer(format(as.Date(events$date), "%m"))
}

periods <- seq(ANALYSIS_YR_LO, ANALYSIS_YR_HI, by = PERIOD_LENGTH)
n_periods <- length(periods)
events <- events |>
  mutate(
    period = cut(year, breaks = c(periods - 0.5, ANALYSIS_YR_HI + 0.5),
                  labels = paste0("P", seq_len(n_periods)),
                  right = TRUE)
  ) |>
  filter(!is.na(period))

# ── 2. 分析 A：时间偏差 — 各 period 繁殖季/非繁殖季数据量 ──────────
month_analysis <- events |>
  group_by(period, month) |>
  summarise(n_records = n(), n_species = n_distinct(species), .groups = "drop") |>
  mutate(
    season = case_when(
      month %in% 4:8  ~ "breeding",
      month %in% c(1:3, 9:12) ~ "non_breeding"
    )
  )

period_summary <- month_analysis |>
  group_by(period, season) |>
  summarise(
    n_records = sum(n_records),
    n_species_month_avg = mean(n_species),
    .groups = "drop"
  ) |>
  pivot_wider(names_from = season, values_from = c(n_records, n_species_month_avg),
              values_fill = 0) |>
  mutate(
    total = n_records_breeding + n_records_non_breeding,
    pct_breeding = 100 * n_records_breeding / total,
    pct_non_breeding = 100 * n_records_non_breeding / total
  )

message("[15b] === Period × Season data volume ===")
for (i in seq_len(nrow(period_summary))) {
  r <- period_summary[i, ]
  message(sprintf("  %s: breeding %d (%.1f%%), non-breeding %d (%.1f%%)",
                  r$period, r$n_records_breeding, r$pct_breeding,
                  r$n_records_non_breeding, r$pct_non_breeding))
}

# ── 3. 分析 B：非繁殖季物种组成差异 ────────────────────────────────
# 哪些物种主要出现在非繁殖季？（候鸟/越冬鸟）
species_season <- events |>
  group_by(species, season = ifelse(month %in% 4:8, "breeding", "non_breeding")) |>
  summarise(n = n(), .groups = "drop") |>
  pivot_wider(names_from = season, values_from = n, values_fill = 0) |>
  mutate(
    total = breeding + non_breeding,
    pct_non_breeding = 100 * non_breeding / total
  ) |>
  arrange(desc(pct_non_breeding))

# 非繁殖季占比 >50% 的物种（越冬/迁徙为主）
non_breeding_dominant <- species_season |> filter(pct_non_breeding > 50)
message(sprintf("[15b] Species with >50%% non-breeding records: %d / %d (%.1f%%)",
                nrow(non_breeding_dominant), nrow(species_season),
                100 * nrow(non_breeding_dominant) / nrow(species_season)))

# 非繁殖季占比 >30% 的物种
moderate_non_breeding <- species_season |> filter(pct_non_breeding > 30)
message(sprintf("[15b] Species with >30%% non-breeding records: %d / %d (%.1f%%)",
                nrow(moderate_non_breeding), nrow(species_season),
                100 * nrow(moderate_non_breeding) / nrow(species_season)))

# ── 4. 分析 C：网格级过滤影响 ───────────────────────────────────────
grid_season <- events |>
  group_by(grid_cell, season = ifelse(month %in% 4:8, "breeding", "non_breeding")) |>
  summarise(n_records = n(), n_species = n_distinct(species), .groups = "drop") |>
  pivot_wider(names_from = season, values_from = c(n_records, n_species),
              values_fill = 0)

# 繁殖季过滤后数据损失严重的网格
grid_loss <- grid_season |>
  mutate(
    total_records = n_records_breeding + n_records_non_breeding,
    pct_lost = 100 * n_records_non_breeding / total_records,
    total_species = pmax(n_species_breeding, n_species_non_breeding),
    species_lost = n_species_non_breeding - n_species_breeding
  ) |>
  arrange(desc(pct_lost))

# 统计损失分布
loss_stats <- grid_loss |>
  summarise(
    median_pct_lost = median(pct_lost, na.rm = TRUE),
    mean_pct_lost = mean(pct_lost, na.rm = TRUE),
    p75_pct_lost = quantile(pct_lost, 0.75, na.rm = TRUE),
    n_grid_over50 = sum(pct_lost > 50),
    n_grid_over30 = sum(pct_lost > 30),
    n_total = n()
  )

message("[15b] === Grid-level data loss from breeding filter ===")
message(sprintf("  Median loss: %.1f%%, Mean: %.1f%%, P75: %.1f%%",
                loss_stats$median_pct_lost, loss_stats$mean_pct_lost,
                loss_stats$p75_pct_lost))
message(sprintf("  Grids losing >50%% data: %d/%d, >30%%: %d/%d",
                loss_stats$n_grid_over50, loss_stats$n_total,
                loss_stats$n_grid_over30, loss_stats$n_total))

# ── 5. 分析 D：繁殖季 vs 全年 — 朴素丰富度趋势对比 ────────────────
# 全年数据的网格丰富度趋势
naive_all_year <- events |>
  group_by(grid_cell, period) |>
  summarise(n_species = n_distinct(species), .groups = "drop") |>
  group_by(grid_cell) |>
  arrange(period) |>
  filter(n() >= 3) |>
  summarise(trend_all = coef(lm(n_species ~ seq_len(n())))[2], .groups = "drop")

# 仅繁殖季
naive_breeding <- events |>
  filter(month %in% BREEDING_MONTHS) |>
  group_by(grid_cell, period) |>
  summarise(n_species = n_distinct(species), .groups = "drop") |>
  group_by(grid_cell) |>
  arrange(period) |>
  filter(n() >= 3) |>
  summarise(trend_breeding = coef(lm(n_species ~ seq_len(n())))[2], .groups = "drop")

trend_comparison <- naive_all_year |>
  inner_join(naive_breeding, by = "grid_cell") |>
  mutate(
    sign_all = sign(trend_all),
    sign_breeding = sign(trend_breeding),
    direction_consistent = sign_all == sign_breeding,
    trend_diff = trend_breeding - trend_all
  )

n_consistent <- sum(trend_comparison$direction_consistent, na.rm = TRUE)
n_total <- nrow(trend_comparison)
message(sprintf("[15b] === Trend direction: breeding vs all-year ==="))
message(sprintf("  Consistent: %d/%d (%.1f%%)",
                n_consistent, n_total, 100 * n_consistent / n_total))
message(sprintf("  Median trend difference: %.3f", median(trend_comparison$trend_diff, na.rm = TRUE)))

# 翻转的网格：方向不一致
flipped <- trend_comparison |> filter(!direction_consistent)
if (nrow(flipped) > 0) {
  message(sprintf("  Direction flipped: %d grids", nrow(flipped)))
  # 翻转网格的空间分布（哪些区域受影响最大）
  flipped_summary <- flipped |>
    mutate(direction = case_when(
      sign_all > 0 & sign_breeding < 0 ~ "all_up_breed_down",
      sign_all < 0 & sign_breeding > 0 ~ "all_down_breed_up",
      TRUE ~ "other"
    )) |>
    count(direction)
  for (i in seq_len(nrow(flipped_summary))) {
    message(sprintf("    %s: %d", flipped_summary$direction[i], flipped_summary$n[i]))
  }
}

# ── 6. 分析 E：按 period 的月份分布变化 ──────────────────────────────
# 早期（P1）vs 晚期（P5）的月份分布是否有系统性偏移
monthly_by_period <- events |>
  group_by(period, month) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(period) |>
  mutate(pct = 100 * n / sum(n)) |>
  ungroup()

# 早期 vs 晚期繁殖季占比
breeding_pct_by_period <- monthly_by_period |>
  filter(month %in% 4:8) |>
  group_by(period) |>
  summarise(pct_breeding = sum(pct), .groups = "drop")

message("[15b] === Breeding season proportion by period ===")
for (i in seq_len(nrow(breeding_pct_by_period))) {
  r <- breeding_pct_by_period[i, ]
  message(sprintf("  %s: %.1f%%", r$period, r$pct_breeding))
}

# ── 7. 输出 ──────────────────────────────────────────────────────────
results <- list(
  period_summary = period_summary,
  species_season_top = head(species_season, 100),
  grid_loss_summary = loss_stats,
  trend_comparison = trend_comparison,
  monthly_by_period = monthly_by_period,
  breeding_pct_by_period = breeding_pct_by_period
)

for (nm in names(results)) {
  write_csv(results[[nm]], v3_file("results", paste0("breeding_sensitivity_", nm)))
}

# ── 8. 决策建议 ──────────────────────────────────────────────────────
message("\n[15b] ========== DECISION RECOMMENDATION ==========")

if (n_consistent / n_total > 0.9) {
  message("→ 趋势方向一致性 >90%：繁殖季过滤对结论影响很小")
  message("→ 建议：保持繁殖季过滤（生态学上更合理）")
} else if (n_consistent / n_total > 0.75) {
  message("→ 趋势方向一致性 75-90%：繁殖季过滤有一定影响")
  message("→ 建议：保持繁殖季过滤，但在论文中报告全年敏感性分析结果")
} else {
  message("→ 趋势方向一致性 <75%：繁殖季过滤显著改变结论")
  message("→ 建议：考虑改用全年数据 + month 检测协变量")
}

if (loss_stats$median_pct_lost > 40) {
  message("→ 中位数数据损失 >40%：繁殖季过滤代价较高")
  message("→ 建议：至少对受影响最大的区域做全年数据分析")
} else {
  message("→ 中位数数据损失 <40%：繁殖季过滤数据代价可接受")
}

message("[15b] ===============================================\n")

log_time("15b", "DONE")
