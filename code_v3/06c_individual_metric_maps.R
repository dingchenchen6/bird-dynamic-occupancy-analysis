#!/usr/bin/env Rscript
## 06c_individual_metric_maps.R  —  每个多样性指标 × 每个时期的独立空间地图
##
## 输入: table_community_metrics_with_cri_*.csv
## 输出: figures_v3/maps_per_metric/fig_map_{metric}_{period}.pdf + .png

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
log_time("06c", "Generating individual metric maps")

# ── 加载数据 ──────────────────────────────────────────────────────────
metric_df <- read_csv_safe(v3_file("results", paste0("table_community_metrics_with_cri_", run_label_extended)))
if (is.null(metric_df)) stop("Community metrics not found.")

grid_sf <- safe_read(v3_file("derived", paste0("china_grid_", GRID_SIZE_KM, "km_v3"), "rds"))
if (is.null(grid_sf)) {
  grid_sf <- safe_read(v3_file("derived", "china_grid_100km_v3", "rds"))
}
china_layers <- china_province_basemap()

# 指标标签
METRIC_LABELS <- c(
  corrected_richness   = "Taxonomic richness",
  shannon              = "Shannon diversity",
  inv_simpson          = "Inverse Simpson",
  trait_volume         = "Functional trait volume",
  trait_dispersion     = "Trait dispersion (FDis)",
  rao_q                = "Rao's Q",
  feve                 = "Functional evenness (FEve)",
  fdiv                 = "Functional divergence (FDiv)",
  fric_prob            = "Functional richness (FRic)",
  fdiv_fund            = "FDiv (fundiversity)",
  feve_fund            = "FEve (fundiversity)",
  fdis_prob            = "FDis (fundiversity)",
  fmpd_prob            = "FMPD",
  raoq_fund            = "RaoQ (fundiversity)",
  cwm_pc1              = "CWM PC1",
  cwm_pc2              = "CWM PC2",
  pd_prob              = "Faith's PD",
  mpd_prob             = "MPD",
  pd_prob_mctavish     = "Faith's PD (McTavish)",
  mpd_prob_mctavish    = "MPD (McTavish)"
)

# 创建输出目录
maps_dir <- file.path(DIRS$figures, "maps_per_metric")
if (!dir.exists(maps_dir)) dir.create(maps_dir, recursive = TRUE)

# 确保 grid_sf 是 sf 对象
if (!inherits(grid_sf, "sf")) {
  grid_sf <- sf::st_as_sf(grid_sf, crs = sf::st_crs(4326))
}

# ── 逐指标 × 逐时期生成地图 ──────────────────────────────────────────
all_metrics <- unique(metric_df$metric)
all_periods <- unique(metric_df$period)

message(sprintf("[06c] Generating %d metrics × %d periods = %d maps",
                length(all_metrics), length(all_periods),
                length(all_metrics) * length(all_periods)))

for (m in all_metrics) {
  m_label <- METRIC_LABELS[m] %||% m
  m_data <- metric_df |>
    filter(metric == m) |>
    inner_join(grid_sf, by = "grid_cell") |>
    sf::st_as_sf()  # 保证连接后是 sf 对象

  if (nrow(m_data) == 0) next

  # ── 单时期地图 ──────────────────────────────────────────────────────
  for (p in all_periods) {
    p_data <- m_data |>
      filter(period == p)
    if (nrow(p_data) == 0) next

    p_safe <- gsub(" ", "_", p)
    fname_base <- paste0("fig_map_", m, "_", p_safe)

    # 使用每个指标自身的 scale（非 z-score），便于跨时期比较
    p_plot <- ggplot(p_data) +
      geom_sf(aes(fill = value_mean),
              colour = alpha("white", 0.1), linewidth = 0.03) +
      geom_sf(data = china_layers, fill = NA, colour = "grey40",
              linewidth = 0.2) +
      scale_fill_viridis_c(name = m_label, option = "D", na.value = "grey90") +
      labs(title = paste(m_label, "—", p),
           subtitle = sprintf("Mean = %.2f, SD = %.2f",
                              mean(p_data$value_mean, na.rm = TRUE),
                              sd(p_data$value_mean, na.rm = TRUE))) +
      coord_sf() +
      theme_minimal(base_size = 9) +
      theme(
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey40"),
        legend.position = "right",
        legend.key.width = unit(0.4, "cm"),
        legend.key.height = unit(1.2, "cm"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
      )

    # 保存 PDF
    ggsave(file.path(maps_dir, paste0(fname_base, ".pdf")),
           p_plot, width = 6, height = 5, dpi = 300)
    # 保存 PNG
    ggsave(file.path(maps_dir, paste0(fname_base, ".png")),
           p_plot, width = 6, height = 5, dpi = 300)
  }

  # ── 多 panel 拼合地图（5个时期一排）─────────────────────────────────
  p_all <- ggplot(m_data) +
    geom_sf(aes(fill = value_mean),
            colour = alpha("white", 0.1), linewidth = 0.03) +
    geom_sf(data = china_layers, fill = NA, colour = "grey40",
            linewidth = 0.2) +
    scale_fill_viridis_c(name = m_label, option = "D", na.value = "grey90") +
    facet_wrap(~ period, ncol = 5) +
    coord_sf() +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.key.width = unit(1.5, "cm"),
      legend.key.height = unit(0.4, "cm"),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      strip.text = element_text(size = 10, face = "bold")
    )

  ggsave(file.path(maps_dir, paste0("fig_map_", m, "_all_periods.pdf")),
         p_all, width = 18, height = 5, dpi = 300)
  ggsave(file.path(maps_dir, paste0("fig_map_", m, "_all_periods.png")),
         p_all, width = 18, height = 5, dpi = 300)

  message(sprintf("  [06c] Done: %s (%d periods + multi-panel)", m, length(all_periods)))
}

# ── 同时生成 temporal dynamic 指标地图 ────────────────────────────────
temporal_df <- read_csv_safe(v3_file("results", paste0("table_temporal_dynamics_summary_", run_label_extended)))
if (!is.null(temporal_df)) {
  temporal_metrics <- unique(temporal_df$metric)
  temporal_labels <- c(
    synchrony = "Synchrony", variance_ratio = "Variance ratio",
    turnover_total = "Turnover (total)", turnover_gain = "Turnover (gain)",
    turnover_loss = "Turnover (loss)"
  )
  t_data <- temporal_df |>
    inner_join(grid_sf, by = "grid_cell") |>
    sf::st_as_sf()

  for (m in temporal_metrics) {
    m_label <- temporal_labels[m] %||% m
    m_sub <- t_data |>
      filter(metric == m)
    if (nrow(m_sub) == 0) next

    fname_base <- paste0("fig_map_temporal_", m)
    p_plot <- ggplot(m_sub) +
      geom_sf(aes(fill = mean),
              colour = alpha("white", 0.1), linewidth = 0.03) +
      geom_sf(data = china_layers, fill = NA, colour = "grey40", linewidth = 0.2) +
      scale_fill_viridis_c(name = m_label, option = "C", na.value = "grey90") +
      labs(title = m_label) +
      coord_sf() +
      theme_minimal(base_size = 9) +
      theme(legend.position = "right", panel.grid = element_blank(),
            axis.text = element_blank(), axis.title = element_blank())

    ggsave(file.path(maps_dir, paste0(fname_base, ".pdf")),
           p_plot, width = 6, height = 5, dpi = 300)
    ggsave(file.path(maps_dir, paste0(fname_base, ".png")),
           p_plot, width = 6, height = 5, dpi = 300)
  }
  message("  [06c] Temporal dynamic maps done.")
}

# ── GIF 动画地图 ──────────────────────────────────────────────────────
# 为每个指标生成时期演变的 GIF 动画
# 如果 gganimate 不可用，则生成辅助脚本供本地运行

has_gganimate <- requireNamespace("gganimate", quietly = TRUE)

if (has_gganimate) {
  message("[06c] Generating animated GIFs with gganimate ...")
  library(gganimate)

  for (m in all_metrics) {
    m_label <- METRIC_LABELS[m] %||% m
    m_data <- metric_df |>
      filter(metric == m) |>
      inner_join(grid_sf, by = "grid_cell") |>
      sf::st_as_sf()
    if (nrow(m_data) == 0) next

    # 确保 period 是有序因子，动画按正确顺序播放
    m_data$period <- factor(m_data$period, levels = all_periods)

    p_anim <- ggplot(m_data) +
      geom_sf(aes(fill = value_mean),
              colour = alpha("white", 0.1), linewidth = 0.03) +
      geom_sf(data = china_layers, fill = NA, colour = "grey40",
              linewidth = 0.2) +
      scale_fill_viridis_c(name = m_label, option = "D", na.value = "grey90") +
      labs(title = paste(m_label, "— {closest_state}")) +
      coord_sf() +
      theme_minimal(base_size = 9) +
      theme(
        plot.title = element_text(size = 11, face = "bold"),
        legend.position = "right",
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(1.2, "cm"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
      ) +
      transition_states(period, transition_length = 2, state_length = 3) +
      ease_aes("linear")

    anim_file <- file.path(maps_dir, paste0("fig_map_", m, "_animated.gif"))
    animate(p_anim, nframes = 60, fps = 8, width = 600, height = 500,
            renderer = gifski_renderer(anim_file))
    message(sprintf("  [06c] GIF saved: %s", basename(anim_file)))
  }
  message("[06c] All animated GIFs done.")
} else {
  message("[06c] gganimate not available — skipping GIF generation on server.")
  message("[06c] To generate GIFs locally, install gganimate + gifski in your local R, then loop over fig_map_{metric}_all_periods.png frames or rerun this script there.")
}

log_time("06c", "DONE: individual metric maps")
message(sprintf("[06c] All maps saved to %s/", maps_dir))
