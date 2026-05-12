#!/usr/bin/env Rscript
## 06_figures_publication.R  —  v3 出版级图集
##
## 将 stage-4/5 输出的所有结果整合为 Nature/Science 出版级图集。
## 全部地图遵守 v3 规则：coord_sf、Nature 主题、灰度+accent 配色、
## Arial 7.5pt、89/120/183mm 宽度。
##
## 输出（figures_v3/）：
##   A. 多维多样性时间切片地图（CRI 带宽副图）
##   B. 群落动态趋势 4-panel + CRI 不确定度地图
##   C. Baselga temporal beta 分解双面板
##   D. Naive vs occupancy-corrected richness 对比
##   E. 多样性时序轨迹（带 95% CRI 阴影）
##   F. 性状 CWM 空间模式（含 diet_specialization + habitat_breadth）

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(forcats); library(ggplot2)
  library(patchwork); library(sf)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_mapping.R"))
source(file.path(CODE_V3, "utils_spatial.R"))

P <- ensure_v3_dirs()
RUN_LABEL <- Sys.getenv("V3_RUN_LABEL", RUN_LABEL)

log_time("06", sprintf("Assembling publication figures for %s", RUN_LABEL))

# ── 共享数据 ──────────────────────────────────────────────────────────
grid_sf <- safe_read(v3_file("derived", paste0("china_grid_", GRID_SIZE_KM, "km_v3"), "rds"))
if (is.null(grid_sf)) {
  grid_sf <- safe_read(v3_file("derived", "china_grid_100km_v3", "rds"))
}
if (is.null(grid_sf)) {
  grid_sf <- safe_read(file.path(DIRS$v2_derived, "china_grid_100km_v2.rds"))
}
primary_blocks <- read_csv_safe(v3_file("results", "table_primary_5year_blocks"))
china_layers <- china_province_basemap()

# ── 指标标签映射（v3 含新性状维度）──────────────────────────────────────
METRIC_LABELS <- c(
  corrected_richness   = "Taxonomic richness",
  shannon              = "Shannon diversity",
  pd_prob              = "Faith's PD (prob-weighted)",
  mpd_prob             = "MPD (prob-weighted)",
  trait_volume         = "Functional trait volume",
  rao_q                = "Rao's Q",
  feve                 = "Functional evenness (FEve)",
  fdiv                 = "Functional divergence (FDiv)",
  cwm_diet_specialization = "CWM diet specialization",
  cwm_habitat_breadth     = "CWM habitat breadth"
)

# 用于4面板趋势图的核心指标
PANEL_METRICS <- c("corrected_richness", "shannon", "trait_volume", "pd_prob")

# ── A. 时间切片多样性地图 ─────────────────────────────────────────────

metrics_path <- v3_file("results", paste0("table_community_metrics_with_cri_", RUN_LABEL))
if (file.exists(metrics_path)) {
  m_long <- read_csv(metrics_path, show_col_types = FALSE)

  # 显示的指标（含新性状 CWM）
  show_metrics <- c("corrected_richness", "shannon", "pd_prob", "mpd_prob",
                    "trait_volume", "rao_q", "feve", "fdiv")
  m_show <- m_long |>
    filter(metric %in% show_metrics) |>
    mutate(
      metric_label = recode(metric, !!!METRIC_LABELS),
      metric_label = factor(metric_label, levels = METRIC_LABELS[show_metrics]),
      block_label  = factor(block_label,
                            levels = primary_blocks$block_label)
    ) |>
    filter(!is.na(metric_label), !is.na(block_label)) |>
    group_by(metric_label) |>
    mutate(value_z = as.numeric(scale(value_mean)),
           value_z = pmin(pmax(value_z, -2.5), 2.5)) |>
    ungroup()
  m_sf <- grid_sf |> inner_join(m_show, by = "grid_cell")

  ts_plot <- ggplot(m_sf) +
    geom_sf(aes(fill = value_z),
            colour = scales::alpha("white", 0.08), linewidth = 0.04) +
    geom_sf(data = china_layers, fill = NA, colour = "grey50",
            linewidth = NATURE_LINE_MAP) +
    scale_fill_nature_c(name = "z-score") +
    facet_grid(metric_label ~ block_label) +
    coord_sf() +
    theme_nature_map() +
    theme(strip.text = element_text(size = 6))
  save_nature(ts_plot,
              paste0("fig_v3_multidiversity_timeslices_", RUN_LABEL),
              width_mm = NATURE_WIDTH_L)

  # ── A'. 不确定度（CRI 宽度）地图 ─────────────────────────────────────
  m_sd <- m_show |>
    mutate(cri_width = value_u95 - value_l95) |>
    group_by(metric_label) |>
    mutate(width_z = as.numeric(scale(cri_width)),
           width_z = pmin(pmax(width_z, -2.5), 2.5)) |>
    ungroup()
  m_sd_sf <- grid_sf |> inner_join(m_sd, by = "grid_cell")

  uncertainty_plot <- ggplot(m_sd_sf) +
    geom_sf(aes(fill = cri_width),
            colour = scales::alpha("white", 0.08), linewidth = 0.04) +
    geom_sf(data = china_layers, fill = NA, colour = "grey50",
            linewidth = NATURE_LINE_MAP) +
    scale_fill_gradient(low = "white", high = NATURE_ACCENT,
                        name = "95% CRI width") +
    facet_grid(metric_label ~ block_label, scales = "fixed") +
    coord_sf() +
    theme_nature_map() +
    theme(strip.text = element_text(size = 6))
  save_nature(uncertainty_plot,
              paste0("fig_v3_multidiversity_uncertainty_", RUN_LABEL),
              width_mm = NATURE_WIDTH_L)
}

# ── B. 群落动态趋势 4-panel + 不确定度 ────────────────────────────────
# FIX #4: 读取 05 实际输出的 table_trend_summary_* 文件
#         并从中构建 sig95 列（95% CRI 不含 0）

trend_path <- v3_file("results", paste0("table_trend_summary_", RUN_LABEL))
if (!file.exists(trend_path)) {
  # fallback: 尝试旧文件名
  trend_path <- v3_file("results", paste0("table_grid_trends_with_cri_", RUN_LABEL))
}
if (file.exists(trend_path)) {
  trd <- read_csv(trend_path, show_col_types = FALSE)

  # 优先使用 theil_sen 方法（更稳健）
  if ("method" %in% names(trd)) {
    trd <- trd |> filter(method == "theil_sen")
  }

  # FIX #4: 如果没有 sig95 列，从 CRI 推断
  if (!"sig95" %in% names(trd)) {
    if (all(c("q025", "q975") %in% names(trd))) {
      trd$sig95 <- (trd$q025 > 0) | (trd$q975 < 0)
    } else {
      trd$sig95 <- FALSE
    }
  }

  # 趋势均值列：可能叫 mean 或 trend_mean
  if (!"trend_mean" %in% names(trd) && "mean" %in% names(trd)) {
    trd$trend_mean <- trd$mean
  }

  trd_show <- trd |>
    filter(metric %in% PANEL_METRICS) |>
    mutate(metric_label = recode(metric,
      corrected_richness = "Richness trend",
      shannon            = "Shannon trend",
      trait_volume       = "Trait-volume trend",
      pd_prob            = "Faith's PD trend")) |>
    group_by(metric_label) |>
    mutate(z = as.numeric(scale(trend_mean)),
           z = pmin(pmax(z, -2.5), 2.5)) |>
    ungroup()
  trd_sf <- grid_sf |> inner_join(trd_show, by = "grid_cell")

  trd_plot <- ggplot(trd_sf) +
    geom_sf(aes(fill = z),
            colour = scales::alpha("white", 0.08), linewidth = 0.04) +
    geom_sf(data = trd_sf |> filter(sig95),
            aes(geometry = geometry),
            inherit.aes = FALSE, fill = NA,
            colour = "#222222", linewidth = 0.10) +
    geom_sf(data = china_layers, fill = NA, colour = "grey50",
            linewidth = NATURE_LINE_MAP) +
    scale_fill_nature_c(name = "Trend z-score") +
    facet_wrap(~ metric_label, ncol = 2) +
    coord_sf() +
    theme_nature_map()
  save_nature(trd_plot,
              paste0("fig_v3_community_trends_with_cri_", RUN_LABEL),
              width_mm = NATURE_WIDTH_L, height_mm = 140)
}

# ── C. Temporal beta：Baselga 双面板地图 ───────────────────────────────
# FIX #5: 使用 05 实际输出的 table_baselga_summary_* 文件

beta_path <- v3_file("results", paste0("table_baselga_summary_", RUN_LABEL))
if (!file.exists(beta_path)) {
  beta_path <- v3_file("results", paste0("table_temporal_beta_with_cri_", RUN_LABEL))
}
if (file.exists(beta_path)) {
  bt <- read_csv(beta_path, show_col_types = FALSE)

  # 将 05 输出的 metric 名映射到显示标签
  # 05 输出: beta_sim → turnover, beta_sne → nestedness
  bt_grid <- bt |>
    filter(metric %in% c("beta_sim", "beta_sne", "beta_sor")) |>
    mutate(metric = recode(metric,
                           beta_sim = "turnover",
                           beta_sne = "nestedness",
                           beta_sor = "bray")) |>
    group_by(grid_cell, metric) |>
    summarise(value = mean(mean, na.rm = TRUE), .groups = "drop") |>
    mutate(metric_label = recode(metric,
      turnover   = "Baselga turnover",
      nestedness = "Baselga nestedness",
      bray       = "Bray-Curtis"))
  bt_sf <- grid_sf |> inner_join(bt_grid, by = "grid_cell")

  bt_plot <- ggplot(bt_sf) +
    geom_sf(aes(fill = value),
            colour = scales::alpha("white", 0.08), linewidth = 0.04) +
    geom_sf(data = china_layers, fill = NA, colour = "grey50",
            linewidth = NATURE_LINE_MAP) +
    scale_fill_gradient(low = "white", high = NATURE_ACCENT,
                        name = "Mean dissimilarity") +
    facet_wrap(~ metric_label, ncol = 3) +
    coord_sf() +
    theme_nature_map()
  save_nature(bt_plot,
              paste0("fig_v3_temporal_beta_baselga_", RUN_LABEL),
              width_mm = NATURE_WIDTH_L, height_mm = 60)
}

# ── D. Naive vs occupancy-corrected richness ───────────────────────────

naive_grid <- read_csv_safe(v3_file("results", "table_species_detection_coverage"))
sv_path <- v3_file("derived", "species_visit_2000_2024", "rds")
if (!is.null(naive_grid) && file.exists(sv_path)) {
  sv <- readRDS(sv_path)
  naive_block <- as_tibble(sv) |>
    group_by(grid_cell, block_id) |>
    summarise(naive_richness = n_distinct(species), .groups = "drop")
  # FIX #11: 用 block_label 替代 block_id 做 join，避免编码不一致
  # primary_blocks 包含 block_id → block_label 映射
  if ("block_label" %in% names(primary_blocks)) {
    naive_block <- naive_block |>
      left_join(primary_blocks |> select(block_id, block_label),
                by = "block_id")
  } else {
    # fallback: block_id 与 block_label 相同
    naive_block$block_label <- naive_block$block_id
  }

  if (file.exists(metrics_path)) {
    corr_block <- m_long |>
      filter(metric == "corrected_richness") |>
      select(grid_cell, block_label, corrected_richness = value_mean)
    cmp <- naive_block |>
      select(grid_cell, block_label, naive_richness) |>
      inner_join(corr_block, by = c("grid_cell", "block_label"))
    cmp_plot <- ggplot(cmp, aes(naive_richness, corrected_richness,
                                 colour = block_label)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2,
                  colour = "grey50", linewidth = 0.2) +
      geom_point(alpha = 0.55, size = 1.0) +
      scale_colour_grey(name = "Period") +
      facet_wrap(~ block_label, ncol = 5) +
      labs(x = "Naive species richness (observed)",
           y = "Occupancy-corrected richness") +
      theme_nature_pub() +
      theme(legend.position = "none")
    save_nature(cmp_plot,
                paste0("fig_v3_naive_vs_corrected_richness_", RUN_LABEL),
                width_mm = NATURE_WIDTH_L, height_mm = 50)
  }
}

# ── E. 多样性时序轨迹（CRI 阴影）──────────────────────────────────────

if (file.exists(metrics_path)) {
  traj_metrics <- c("corrected_richness", "shannon",
                    "pd_prob", "trait_volume", "rao_q",
                    "feve", "fdiv",
                    "cwm_diet_specialization", "cwm_habitat_breadth")
  traj <- m_long |>
    filter(metric %in% traj_metrics) |>
    group_by(metric, block_label) |>
    summarise(
      mean_med = median(value_mean, na.rm = TRUE),
      l95 = quantile(value_mean, 0.025, na.rm = TRUE),
      u95 = quantile(value_mean, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      block_label = factor(block_label,
                           levels = primary_blocks$block_label),
      metric_label = recode(metric, !!!METRIC_LABELS)
    ) |>
    filter(!is.na(metric_label))

  traj_plot <- ggplot(traj, aes(block_label, mean_med,
                                 group = metric_label)) +
    geom_ribbon(aes(ymin = l95, ymax = u95),
                alpha = 0.15, fill = NATURE_ACCENT) +
    geom_line(colour = NATURE_ACCENT, linewidth = 0.4) +
    geom_point(colour = NATURE_ACCENT, size = 1.2) +
    facet_wrap(~ metric_label, scales = "free_y", ncol = 3) +
    labs(x = NULL, y = NULL) +
    theme_nature_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_nature(traj_plot,
              paste0("fig_v3_multidiversity_trajectories_", RUN_LABEL),
              width_mm = NATURE_WIDTH_L, height_mm = 100)
}

# ── F. 性状 CWM 空间模式（v3 新增：含 diet_specialization + habitat_breadth）──

trait_cwm_path <- v3_file("results",
                           paste0("table_cwm_spatial_pattern_", RUN_LABEL))
if (file.exists(trait_cwm_path)) {
  cwm <- read_csv(trait_cwm_path, show_col_types = FALSE)
  cwm_long <- cwm |>
    pivot_longer(matches("^(cwm_|z_)"), names_to = "trait", values_to = "value") |>
    filter(trait %in% c("cwm_diet_specialization", "cwm_habitat_breadth")) |>
    mutate(trait_label = recode(trait,
      cwm_diet_specialization = "CWM diet specialization",
      cwm_habitat_breadth     = "CWM habitat breadth"))
  cwm_sf <- grid_sf |> inner_join(cwm_long, by = "grid_cell")

  if (nrow(cwm_sf) > 0) {
    cwm_plot <- ggplot(cwm_sf) +
      geom_sf(aes(fill = value),
              colour = scales::alpha("white", 0.08), linewidth = 0.04) +
      geom_sf(data = china_layers, fill = NA, colour = "grey50",
              linewidth = NATURE_LINE_MAP) +
      scale_fill_nature_c(name = "CWM value", palette = "accent") +
      facet_wrap(~ trait_label, ncol = 2) +
      coord_sf() +
      theme_nature_map()
    save_nature(cwm_plot,
                paste0("fig_v3_cwm_spatial_traits_", RUN_LABEL),
                width_mm = NATURE_WIDTH_M, height_mm = 70)
  }
}

# ── G. 物种水平趋势系统图（v3 新增 Section 10）────────────────────────

# G1. 扩张/稳定/收缩物种分组柱状图
sp_classify_path <- v3_file("results", paste0("table_species_trend_classify_", RUN_LABEL))
if (file.exists(sp_classify_path)) {
  sp_cls <- read_csv(sp_classify_path, show_col_types = FALSE)

  cls_counts <- sp_cls |>
    count(trend_class, name = "n") |>
    mutate(trend_class = factor(trend_class,
                                levels = c("contracting", "stable", "expanding"),
                                labels = c("Contracting", "Stable", "Expanding")))

  cls_plot <- ggplot(cls_counts, aes(trend_class, n, fill = trend_class)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = n), vjust = -0.5, size = 7.5 / 2.835) +
    scale_fill_manual(values = c("Contracting" = "#B2182B",
                                  "Stable"     = "#969696",
                                  "Expanding"  = "#2166AC"),
                      name = NULL) +
    labs(x = NULL, y = "Number of species") +
    theme_nature_pub() +
    theme(legend.position = "none")
  save_nature(cls_plot,
              paste0("fig_v3_species_trend_classification_", RUN_LABEL),
              width_mm = NATURE_WIDTH_S, height_mm = 70)

  # G2. 物种趋势分布密度图
  sp_trend_path <- v3_file("results", paste0("table_species_trend_", RUN_LABEL))
  if (file.exists(sp_trend_path)) {
    sp_tr <- read_csv(sp_trend_path, show_col_types = FALSE) |>
      filter(method == "theil_sen")

    dens_plot <- ggplot(sp_tr, aes(x = mean)) +
      geom_density(fill = NATURE_ACCENT, alpha = 0.3, colour = NATURE_ACCENT) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey50", linewidth = 0.2) +
      labs(x = "Occupancy trend (Theil-Sen slope)", y = "Density") +
      theme_nature_pub()
    save_nature(dens_plot,
                paste0("fig_v3_species_trend_density_", RUN_LABEL),
                width_mm = NATURE_WIDTH_S, height_mm = 60)
  }
}

# G3. 扩张/收缩热点地图
hotspot_path <- v3_file("results", paste0("table_species_hotspot_", RUN_LABEL))
if (file.exists(hotspot_path)) {
  hs <- read_csv(hotspot_path, show_col_types = FALSE)
  hs_sf <- grid_sf |> inner_join(hs, by = "grid_cell")

  if (nrow(hs_sf) > 0) {
    # FIX #10: 先过滤 NA，再安全提取最终 period
    hs_period_ok <- hs_sf |>
      filter(!is.na(period))
    if (nrow(hs_period_ok) > 0) {
      max_p <- max(as.numeric(gsub("P", "", hs_period_ok$period)), na.rm = TRUE)
      hs_final <- hs_period_ok |> filter(period == paste0("P", max_p))
    } else {
      hs_final <- hs_period_ok
    }

    if (nrow(hs_final) > 0) {
      hotspot_plot <- ggplot(hs_final) +
        geom_sf(aes(fill = net_trend_index),
                colour = scales::alpha("white", 0.08), linewidth = 0.04) +
        geom_sf(data = china_layers, fill = NA, colour = "grey50",
                linewidth = NATURE_LINE_MAP) +
        scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC",
                             midpoint = 0, name = "Net trend index") +
        coord_sf() +
        theme_nature_map()
      save_nature(hotspot_plot,
                  paste0("fig_v3_species_hotspot_map_", RUN_LABEL),
                  width_mm = NATURE_WIDTH_M, height_mm = 80)
    }
  }
}

# ── H. Baselga 分解比例图（turnover vs nestedness）───────────────────

baselga_global_path <- v3_file("results", paste0("table_baselga_global_", RUN_LABEL))
if (file.exists(baselga_global_path)) {
  bg <- read_csv(baselga_global_path, show_col_types = FALSE)

  bg_long <- bg |>
    select(period_pair, prop_turnover_mean) |>
    mutate(prop_nestedness = 1 - prop_turnover_mean) |>
    pivot_longer(cols = c(prop_turnover_mean, prop_nestedness),
                 names_to = "component", values_to = "proportion") |>
    mutate(component = recode(component,
                              prop_turnover_mean = "Turnover",
                              prop_nestedness = "Nestedness"))

  baselga_bar_plot <- ggplot(bg_long, aes(period_pair, proportion, fill = component)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = c("Turnover" = NATURE_ACCENT,
                                  "Nestedness" = "#D4DEE4"),
                      name = NULL) +
    labs(x = "Period pair", y = "Proportion of β diversity") +
    theme_nature_pub() +
    theme(legend.position = "top")
  save_nature(baselga_bar_plot,
              paste0("fig_v3_baselga_proportion_", RUN_LABEL),
              width_mm = NATURE_WIDTH_M, height_mm = 70)
}

log_time("06", "Publication figures assembled.")
