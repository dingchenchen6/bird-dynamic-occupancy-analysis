#!/usr/bin/env Rscript
## 06b_regenerate_driver_plots.R  —  v3 驱动因子 + RF 重要性图
##
## 从 brms driver fit 和 RF 重要性结果生成 Nature 风格图表：
##   1. 变化量驱动 raincloud 图（brms 后验分布，优先）
##   2. 多响应汇总 facet 图
##   3. RF 排列重要性柱状图（带后验 CI，变化量驱动）
##   4. varpart vs RF 对比面板图
##   5. 敏感性分析对比图
##
## v3 改进：优先使用变化量驱动（DRIVER_GROUPS_TREND），
##   静态驱动（DRIVER_GROUPS）作为补充

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(patchwork); library(forcats); library(stringr)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_mapping.R"))
source(file.path(CODE_V3, "utils_importance.R"))
source(file.path(CODE_V3, "utils_plots_advanced.R"))

P <- ensure_v3_dirs()
RUN_LABEL <- Sys.getenv("V3_RUN_LABEL", RUN_LABEL)

log_time("06b", sprintf("Regenerating driver + RF plots for %s (with temporal change drivers)", RUN_LABEL))

# ── 通用标签 ──────────────────────────────────────────────────────────
RESP_PRETTY <- c(
  trend_corrected_richness = "Richness trend",
  trend_shannon            = "Shannon trend",
  trend_pd_prob            = "Faith's PD trend",
  trend_trait_volume       = "Trait-volume trend",
  trend_rao_q              = "Rao's Q trend",
  trend_mpd_prob           = "MPD trend"
)

# ── 1. brms 驱动因子 raincloud 图 ─────────────────────────────────────

# 优先加载变化量驱动的 brms fit
fit_files <- list.files(DIRS$derived,
  pattern = sprintf("^brms_driver_trend_.*_%s\\.rds$", RUN_LABEL),
  full.names = TRUE)

# 如果没有变化量驱动的 fit，尝试静态版
if (length(fit_files) == 0) {
  fit_files <- list.files(DIRS$derived,
    pattern = sprintf("^brms_driver_fit_.*_%s\\.rds$", RUN_LABEL),
    full.names = TRUE)
  driver_type <- "static"
} else {
  driver_type <- "temporal_change"
}

message(sprintf("[06b] Found %d brms fits (type: %s)", length(fit_files), driver_type))

if (length(fit_files) > 0) {
  resp_of <- function(p) {
    bn <- basename(p)
    bn <- sub(paste0("_", RUN_LABEL, "\\.rds$"), "", bn)
    bn <- sub("^brms_driver_trend_", "", bn)
    bn <- sub("^brms_driver_fit_", "", bn)
    bn
  }
  fits <- lapply(fit_files, readRDS)
  names(fits) <- vapply(fit_files, resp_of, character(1))

  title_for <- function(nm) {
    if (nm %in% names(RESP_PRETTY)) RESP_PRETTY[[nm]] else nm
  }

  # 单响应 raincloud
  for (nm in names(fits)) {
    draws <- brms_fixef_draws(fits[[nm]])
    draws$term <- gsub("^z_|^scale", "", draws$term)
    # 应用变化量驱动标签
    draws$term <- recode(draws$term, !!!DRIVER_TREND_LABELS)

    p <- raincloud_nature(draws, x_var = value, y_var = term) +
      geom_vline(xintercept = 0, linetype = "dashed",
                 colour = "grey50", linewidth = 0.2) +
      labs(x = "Standardised coefficient", y = NULL,
           title = if (driver_type == "temporal_change")
             "Drivers of richness trend (temporal change)" else NULL) +
      theme_nature_pub()
    save_nature(p,
                sprintf("fig_v3_driver_raincloud_%s_%s", nm, RUN_LABEL),
                width_mm = NATURE_WIDTH_M, height_mm = 80)
  }
  log_time("06b", sprintf("%d single-response raincloud figures done", length(fits)))

  # 多响应汇总
  draws_all <- purrr::imap_dfr(fits, function(fit, nm) {
    d <- brms_fixef_draws(fit)
    d$term <- gsub("^z_|^scale", "", d$term)
    d$term <- recode(d$term, !!!DRIVER_TREND_LABELS)
    d$response <- title_for(nm)
    d
  })
  med_within <- draws_all |>
    group_by(response, term) |>
    summarise(med = median(value), .groups = "drop")
  draws_all <- draws_all |>
    left_join(med_within, by = c("response", "term")) |>
    mutate(term = forcats::fct_reorder(term, med))

  multi_plot <- ggplot(draws_all, aes(value, term)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey55", linewidth = 0.2) +
    ggdist::stat_halfeye(.width = c(0.5, 0.95), thickness = 0.55,
                         slab_alpha = 0.55, slab_colour = NA,
                         interval_colour = "grey25", point_size = 1.2,
                         fill = NATURE_ACCENT) +
    facet_wrap(~ response, scales = "free", ncol = 2) +
    labs(x = "Standardised coefficient", y = NULL) +
    theme_nature_pub() +
    theme(panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.15))
  save_nature(multi_plot,
              sprintf("fig_v3_driver_raincloud_multipanel_%s", RUN_LABEL),
              width_mm = NATURE_WIDTH_L, height_mm = 140)
  log_time("06b", "Multi-response raincloud figure done")

} else {
  log_time("06b", "No brms driver fits found, skipping raincloud plots")
}

# ── 2. RF 排列重要性图（变化量驱动版）─────────────────────────────────

# 优先使用变化量驱动的 RF 结果
rf_summary_path <- v3_file("results", "table_rf_importance_variable_summary_trend")
rf_group_path   <- v3_file("results", "table_rf_importance_group_summary_trend")

# 回退到旧版
if (!file.exists(rf_summary_path)) {
  rf_summary_path <- v3_file("results", "table_rf_importance_variable_summary")
  rf_group_path   <- v3_file("results", "table_rf_importance_group_summary")
}

if (file.exists(rf_summary_path)) {
  rf_summary <- read_csv(rf_summary_path, show_col_types = FALSE)

  # 变量重要性柱状图（带 95% CI）
  rf_var_plot <- plot_rf_variable_importance(rf_summary, top_n = 15) +
    labs(x = NULL, y = "Permutation importance")
  save_nature(rf_var_plot,
              paste0("fig_v3_rf_variable_importance_", RUN_LABEL),
              width_mm = NATURE_WIDTH_S, height_mm = 80)
  log_time("06b", "RF variable importance figure done")

  # 组级重要性柱状图
  if (file.exists(rf_group_path)) {
    rf_group <- read_csv(rf_group_path, show_col_types = FALSE)

    # 应用组标签
    if ("group" %in% names(rf_group)) {
      rf_group <- rf_group |>
        mutate(group_label = recode(group, !!!DRIVER_GROUP_TREND_LABELS))
    }

    rf_grp_plot <- ggplot(rf_group,
                           aes(x = reorder(
                             if ("group_label" %in% names(rf_group)) group_label else group,
                             mean), y = mean)) +
      geom_col(fill = NATURE_ACCENT, alpha = 0.85, width = 0.55) +
      geom_errorbar(aes(ymin = q_lo, ymax = q_hi), width = 0.2,
                    linewidth = 0.3) +
      coord_flip() +
      labs(x = NULL, y = "Group importance (mean +/- 95% CI)") +
      theme_nature_pub()
    save_nature(rf_grp_plot,
                paste0("fig_v3_rf_group_importance_", RUN_LABEL),
                width_mm = NATURE_WIDTH_S, height_mm = 55)
    log_time("06b", "RF group importance figure done")
  }

  # ── 3. varpart vs RF 对比面板图 ──────────────────────────────────────
  varpart_path <- v3_file("results",
                           paste0("table_varpart_richness_trend_", RUN_LABEL))
  if (file.exists(varpart_path) && file.exists(rf_group_path)) {
    varpart_df <- read_csv(varpart_path, show_col_types = FALSE)
    vp_summary <- varpart_df |>
      filter(grepl("pure", component, ignore.case = TRUE)) |>
      mutate(
        # 使用变化量驱动标签
        component = recode(component,
          `Climate change (pure)`    = "Climate change",
          `Land use change (pure)`   = "Land use change",
          `Human pressure change (pure)` = "Human pressure change",
          `Spatial baseline (pure)`  = "Spatial baseline",
          `Climate (pure)`     = "Climate",
          `Topo+Habitat (pure)` = "Topo+Habitat",
          `Human (pure)`       = "Human",
          `Space (pure)`       = "Space"
        ),
        adj_R2 = as.numeric(adj_R2)
      ) |>
      filter(!is.na(adj_R2))

    compare_plot <- plot_rf_vs_varpart(vp_summary, rf_group)
    save_nature(compare_plot,
                paste0("fig_v3_varpart_vs_rf_importance_", RUN_LABEL),
                width_mm = NATURE_WIDTH_L, height_mm = 60)
    log_time("06b", "varpart vs RF comparison figure done")
  }
} else {
  log_time("06b", "No RF importance results found, skipping RF figures")
}

# ── 4. 敏感性分析对比图 ───────────────────────────────────────────────
sens_path <- v3_file("results",
  paste0("table_varpart_richness_trend_sensitivity_", RUN_LABEL))
if (file.exists(sens_path)) {
  sens_df <- read_csv(sens_path, show_col_types = FALSE)
  log_time("06b", "Sensitivity varpart data available — plot to be generated")
  # TODO: 绘制主模型 vs 敏感性分析的对比面板图
}

log_time("06b", "Done.")
