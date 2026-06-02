#!/usr/bin/env Rscript
## 06d_stability_figures.R  —  群落稳定性可视化
##
## Scientific question / 科学问题:
## What are the spatial patterns of community stability, resilience, and resistance?
## How do they relate to environmental drivers?
## 群落稳定性、韧性和抗性的空间格局是什么？它们与环境驱动因子有何关系？
##
## Objective / 分析目标:
## Generate publication-ready maps and scatter plots for stability metrics.
## 生成可用于论文的稳定性指标地图和散点图。
##
## Input / 输入:
##   - table_stability_summary_*.csv
##   - table_resistance_summary_*.csv
##   - grid_environment_*.rds (for driver variables)
##   - china_grid_*.rds (spatial grid)
##
## Output / 输出:
##   - figures_v3/fig_stability_map_*.pdf/png
##   - figures_v3/fig_stability_vs_drivers_*.pdf/png
##   - figures_v3/fig_resilience_resistance_scatter.pdf/png

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(sf); library(scales)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_mapping.R"))
P <- ensure_v3_dirs()

run_label_extended <- paste0(RUN_LABEL, "_extended")
log_time("06d", "Starting stability figures")

# ── 加载数据 ──────────────────────────────────────────────────────────
stability_df <- read_csv_safe(v3_file("results", paste0("table_stability_summary_", run_label_extended)))
resistance_df <- read_csv_safe(v3_file("results", paste0("table_resistance_summary_", run_label_extended)))

grid_sf <- safe_read(v3_file("derived", paste0("china_grid_", GRID_SIZE_KM, "km_v3"), "rds"))
if (is.null(grid_sf)) {
  grid_sf <- safe_read(v3_file("derived", "china_grid_100km_v3", "rds"))
}

china_layers <- china_province_basemap()

# 确保 grid_sf 是 sf 对象
if (!is.null(grid_sf) && !inherits(grid_sf, "sf")) {
  grid_sf <- sf::st_as_sf(grid_sf, crs = sf::st_crs(4326))
}

# ── 指标标签映射 ──────────────────────────────────────────────────────
STABILITY_LABELS <- c(
  temporal_stability_cv_inv  = "Temporal stability (CV⁻¹)",
  temporal_stability_logcv   = "log(CV) — lower = stable",
  trajectory_distance        = "Composition trajectory distance",
  mean_bray_curtis_change    = "Mean Bray-Curtis change",
  ar1_resilience             = "Resilience (AR1 coefficient)"
)

RESISTANCE_LABELS <- c(
  hfi_change     = "Resistance to HFI change",
  climate_change = "Resistance to climate change"
)

METRIC_LABELS_SHORT <- c(
  corrected_richness   = "Richness",
  shannon              = "Shannon",
  inv_simpson          = "Inv. Simpson",
  trait_volume         = "Trait volume",
  trait_dispersion     = "FDis",
  rao_q                = "Rao's Q",
  feve                 = "FEve",
  fdiv                 = "FDiv",
  fric_prob            = "FRic",
  fdiv_fund            = "FDiv (fundiv)",
  feve_fund            = "FEve (fundiv)",
  fdis_prob            = "FDis (fundiv)",
  fmpd_prob            = "FMPD",
  raoq_fund            = "RaoQ (fundiv)",
  cwm_pc1              = "CWM PC1",
  cwm_pc2              = "CWM PC2",
  pd_prob              = "Faith's PD",
  mpd_prob             = "MPD",
  pd_prob_mctavish     = "PD (McTavish)",
  mpd_prob_mctavish    = "MPD (McTavish)"
)

maps_dir <- file.path(DIRS$figures, "stability_maps")
if (!dir.exists(maps_dir)) dir.create(maps_dir, recursive = TRUE)

# ── 1. 稳定性空间地图 ─────────────────────────────────────────────────
if (!is.null(stability_df) && nrow(stability_df) > 0 && !is.null(grid_sf)) {
  message("[06d] Generating stability maps ...")

  # 选择代表性的 metric 组合，避免过多图
  priority_metrics <- c("corrected_richness", "shannon", "trait_volume",
                        "fric_prob", "pd_prob", "mpd_prob")
  available_metrics <- intersect(priority_metrics, unique(stability_df$metric))
  if (length(available_metrics) == 0) {
    available_metrics <- unique(stability_df$metric)[1:min(3, length(unique(stability_df$metric)))]
  }

  stab_metrics <- unique(stability_df$stability_metric)

  for (sm in stab_metrics) {
    sm_label <- STABILITY_LABELS[sm] %||% sm
    sm_data <- stability_df |>
      filter(stability_metric == sm, metric %in% available_metrics) |>
      inner_join(grid_sf, by = "grid_cell") |>
      sf::st_as_sf()

    if (nrow(sm_data) == 0) next

    # facet by metric
    sm_data$metric_lab <- METRIC_LABELS_SHORT[sm_data$metric] %||% sm_data$metric

    p_map <- ggplot(sm_data) +
      geom_sf(aes(fill = mean),
              colour = alpha("white", 0.1), linewidth = 0.03) +
      geom_sf(data = china_layers, fill = NA, colour = "grey40",
              linewidth = 0.2) +
      scale_fill_viridis_c(name = sm_label, option = "D", na.value = "grey90") +
      facet_wrap(~ metric_lab, ncol = 3) +
      coord_sf() +
      theme_minimal(base_size = 9) +
      theme(
        plot.title = element_text(size = 10, face = "bold"),
        legend.position = "bottom",
        legend.key.width = unit(1.2, "cm"),
        legend.key.height = unit(0.3, "cm"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(size = 8, face = "bold")
      )

    fname <- paste0("fig_stability_map_", sm)
    ggsave(file.path(maps_dir, paste0(fname, ".pdf")),
           p_map, width = 10, height = 7, dpi = 300)
    ggsave(file.path(maps_dir, paste0(fname, ".png")),
           p_map, width = 10, height = 7, dpi = 300)
    message(sprintf("  [06d] Saved: %s", fname))
  }
} else {
  message("[06d] Stability data or grid missing — skipping stability maps.")
}

# ── 2. 抗性空间地图 ───────────────────────────────────────────────────
if (!is.null(resistance_df) && nrow(resistance_df) > 0 && !is.null(grid_sf)) {
  message("[06d] Generating resistance maps ...")

  priority_metrics <- c("corrected_richness", "shannon", "trait_volume",
                        "fric_prob", "pd_prob", "mpd_prob")
  available_metrics <- intersect(priority_metrics, unique(resistance_df$metric))
  if (length(available_metrics) == 0) {
    available_metrics <- unique(resistance_df$metric)[1:min(3, length(unique(resistance_df$metric)))]
  }

  pert_types <- unique(resistance_df$perturbation)

  for (pert in pert_types) {
    pert_label <- RESISTANCE_LABELS[pert] %||% pert
    r_data <- resistance_df |>
      filter(perturbation == pert, metric %in% available_metrics) |>
      inner_join(grid_sf, by = "grid_cell") |>
      sf::st_as_sf()

    if (nrow(r_data) == 0) next

    r_data$metric_lab <- METRIC_LABELS_SHORT[r_data$metric] %||% r_data$metric

    p_map <- ggplot(r_data) +
      geom_sf(aes(fill = mean),
              colour = alpha("white", 0.1), linewidth = 0.03) +
      geom_sf(data = china_layers, fill = NA, colour = "grey40",
              linewidth = 0.2) +
      scale_fill_viridis_c(name = pert_label, option = "C", na.value = "grey90") +
      facet_wrap(~ metric_lab, ncol = 3) +
      coord_sf() +
      theme_minimal(base_size = 9) +
      theme(
        plot.title = element_text(size = 10, face = "bold"),
        legend.position = "bottom",
        legend.key.width = unit(1.2, "cm"),
        legend.key.height = unit(0.3, "cm"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(size = 8, face = "bold")
      )

    fname <- paste0("fig_resistance_map_", pert)
    ggsave(file.path(maps_dir, paste0(fname, ".pdf")),
           p_map, width = 10, height = 7, dpi = 300)
    ggsave(file.path(maps_dir, paste0(fname, ".png")),
           p_map, width = 10, height = 7, dpi = 300)
    message(sprintf("  [06d] Saved: %s", fname))
  }
} else {
  message("[06d] Resistance data or grid missing — skipping resistance maps.")
}

# ── 3. 稳定性 vs 驱动因子关系 ─────────────────────────────────────────
grid_env <- safe_read(v3_file("derived", paste0("grid_environment", GRID_TAG, "_v3"), "rds"))
if (is.null(grid_env)) {
  grid_env <- safe_read(v3_file("derived", "grid_environment_v3", "rds"))
}

if (!is.null(stability_df) && !is.null(grid_env) && nrow(stability_df) > 0) {
  message("[06d] Generating stability vs drivers scatter ...")

  # 选择驱动变量
  driver_vars <- c("delta_hfi", "delta_t_mean", "delta_forest",
                   "elev_mean", "bio11")
  avail_drivers <- intersect(driver_vars, names(grid_env))

  if (length(avail_drivers) > 0) {
    # 为每个 stability_metric 和代表性 metric 绘制 scatter
    stab_metrics <- unique(stability_df$stability_metric)
    rep_metric <- "corrected_richness"
    if (!(rep_metric %in% stability_df$metric)) {
      rep_metric <- unique(stability_df$metric)[1]
    }

    for (sm in stab_metrics) {
      sm_data <- stability_df |>
        filter(stability_metric == sm, metric == rep_metric) |>
        inner_join(grid_env, by = "grid_cell")
      if (nrow(sm_data) == 0) next

      sm_label <- STABILITY_LABELS[sm] %||% sm

      # 收集所有 driver 的 long format
      plot_list <- list()
      for (dv in avail_drivers) {
        if (!dv %in% names(sm_data)) next
        dv_label <- DRIVER_TREND_LABELS[dv] %||% dv
        # 去除 NA
        sub <- sm_data[!is.na(sm_data[[dv]]), ]
        if (nrow(sub) < 10) next

        p <- ggplot(sub, aes(x = .data[[dv]], y = mean)) +
          geom_point(alpha = 0.3, size = 0.8, colour = NATURE_ACCENT) +
          geom_smooth(method = "lm", se = TRUE, colour = "firebrick",
                      linewidth = 0.5, fill = "grey80", alpha = 0.3) +
          labs(x = dv_label, y = sm_label) +
          theme_nature_pub()
        plot_list[[dv]] <- p
      }

      if (length(plot_list) > 0) {
        n <- length(plot_list)
        ncol <- min(3, n)
        nrow <- ceiling(n / ncol)
        # 使用 patchwork 或 cowplot 拼接；如无则逐个保存
        if (requireNamespace("patchwork", quietly = TRUE)) {
          library(patchwork)
          p_comb <- wrap_plots(plot_list, ncol = ncol) +
            plot_annotation(title = paste(sm_label, "—", METRIC_LABELS_SHORT[rep_metric]),
                            theme = theme(plot.title = element_text(size = 11, face = "bold")))
          fname <- paste0("fig_stability_vs_drivers_", sm)
          ggsave(file.path(DIRS$figures, paste0(fname, ".pdf")),
                 p_comb, width = 3 * ncol + 1, height = 3 * nrow, dpi = 300)
          ggsave(file.path(DIRS$figures, paste0(fname, ".png")),
                 p_comb, width = 3 * ncol + 1, height = 3 * nrow, dpi = 300)
          message(sprintf("  [06d] Saved: %s", fname))
        } else {
          # 逐个保存
          for (dv in names(plot_list)) {
            fname <- paste0("fig_stability_vs_drivers_", sm, "_", dv)
            ggsave(file.path(DIRS$figures, paste0(fname, ".pdf")),
                   plot_list[[dv]], width = 5, height = 4, dpi = 300)
            ggsave(file.path(DIRS$figures, paste0(fname, ".png")),
                   plot_list[[dv]], width = 5, height = 4, dpi = 300)
          }
          message(sprintf("  [06d] Saved individual driver plots for %s", sm))
        }
      }
    }
  } else {
    message("[06d] No driver variables available in grid_env.")
  }
} else {
  message("[06d] Skipping stability vs drivers (data missing).")
}

# ── 4. 韧性-抗性散点图 ────────────────────────────────────────────────
if (!is.null(stability_df) && !is.null(resistance_df) &&
    nrow(stability_df) > 0 && nrow(resistance_df) > 0) {
  message("[06d] Generating resilience-resistance scatter ...")

  # 提取韧性 (AR1) 和抗性 (hfi_change)
  resilience_data <- stability_df |>
    filter(stability_metric == "ar1_resilience")

  resistance_data <- resistance_df |>
    filter(perturbation == "hfi_change")

  if (nrow(resilience_data) > 0 && nrow(resistance_data) > 0) {
    # 取共同 metric
    common_metrics <- intersect(unique(resilience_data$metric),
                                unique(resistance_data$metric))
    if (length(common_metrics) == 0) {
      common_metrics <- unique(resilience_data$metric)[1:min(3, length(unique(resilience_data$metric)))]
    }

    rr_list <- list()
    for (m in common_metrics) {
      res_df <- resilience_data |>
        filter(metric == m) |>
        select(grid_cell, resilience = mean)
      rst_df <- resistance_data |>
        filter(metric == m) |>
        select(grid_cell, resistance = mean)
      rr <- inner_join(res_df, rst_df, by = "grid_cell")
      if (nrow(rr) < 10) next

      m_lab <- METRIC_LABELS_SHORT[m] %||% m
      p <- ggplot(rr, aes(x = resilience, y = resistance)) +
        geom_point(alpha = 0.3, size = 0.8, colour = NATURE_ACCENT) +
        geom_smooth(method = "lm", se = TRUE, colour = "firebrick",
                    linewidth = 0.5, fill = "grey80", alpha = 0.3) +
        labs(
          x = "Resilience (AR1 coefficient)",
          y = "Resistance (|Δdiversity| under HFI change)",
          title = m_lab
        ) +
        theme_nature_pub()
      rr_list[[m]] <- p
    }

    if (length(rr_list) > 0) {
      n <- length(rr_list)
      ncol <- min(3, n)
      nrow <- ceiling(n / ncol)
      if (requireNamespace("patchwork", quietly = TRUE)) {
        library(patchwork)
        p_comb <- wrap_plots(rr_list, ncol = ncol) +
          plot_annotation(
            title = "Resilience vs Resistance",
            theme = theme(plot.title = element_text(size = 11, face = "bold"))
          )
        ggsave(file.path(DIRS$figures, "fig_resilience_resistance_scatter.pdf"),
               p_comb, width = 3 * ncol + 1, height = 3 * nrow, dpi = 300)
        ggsave(file.path(DIRS$figures, "fig_resilience_resistance_scatter.png"),
               p_comb, width = 3 * ncol + 1, height = 3 * nrow, dpi = 300)
      } else {
        for (m in names(rr_list)) {
          fname <- paste0("fig_resilience_resistance_", m)
          ggsave(file.path(DIRS$figures, paste0(fname, ".pdf")),
                 rr_list[[m]], width = 5, height = 4, dpi = 300)
          ggsave(file.path(DIRS$figures, paste0(fname, ".png")),
                 rr_list[[m]], width = 5, height = 4, dpi = 300)
        }
      }
      message("  [06d] Saved resilience-resistance scatter plot(s).")
    }
  } else {
    message("[06d] Resilience or resistance data empty — skipping scatter.")
  }
} else {
  message("[06d] Skipping resilience-resistance scatter (data missing).")
}

log_time("06d", "DONE: stability figures")
message(sprintf("[06d] All stability figures saved to %s/ and %s/",
                DIRS$figures, maps_dir))
