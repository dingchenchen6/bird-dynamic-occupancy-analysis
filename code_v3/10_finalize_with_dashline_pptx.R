#!/usr/bin/env Rscript
## 10_finalize_with_dashline_pptx.R  —  v3 PPTX 图集
##
## 1) 重画全部核心图（Nature 风格），地图统一加十段线
## 2) 新增 RF 重要性图
## 3) 集成 PPTX deck（地图用 PNG 嵌入、统计图用 vector dml）
##
## 输出：figures_v3/*_v3.{png,pdf,pptx} + master deck

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(stringr); library(forcats); library(ggplot2); library(patchwork); library(sf)
  library(ggridges)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_spatial.R"))
source(file.path(CODE_V3, "utils_mapping.R"))
source(file.path(CODE_V3, "utils_importance.R"))
source(file.path(CODE_V3, "utils_plots_advanced.R"))

P <- ensure_v3_dirs()
RUN_LABEL <- Sys.getenv("V3_RUN_LABEL", RUN_LABEL)

log_time("10", sprintf("Regenerating maps + PPTX for %s", RUN_LABEL))

# ── 共享数据 ──────────────────────────────────────────────────────────
china_layers <- china_province_basemap()
grid_sf <- readRDS(v3_file("derived", "china_grid_100km_v3", "rds"))
primary_blocks <- read_csv_safe(v3_file("results", "table_primary_5year_blocks"))

# 两个 deck：地图（PNG 嵌入）、统计图（vector dml）
deck_maps  <- list()
deck_stats <- list()

# ── 辅助函数 ──────────────────────────────────────────────────────────
add_plot <- function(name, plot, width_mm = NATURE_WIDTH_L, height_mm = NULL,
                     is_map = FALSE) {
  save_nature(plot, paste0(name, "_v3"), width_mm = width_mm,
              height_mm = height_mm)
  if (is_map) {
    deck_maps[[name]] <<- plot
  } else {
    deck_stats[[name]] <<- plot
    # 单独可编辑 PPTX
    if (requireNamespace("officer", quietly = TRUE) &&
        requireNamespace("rvg", quietly = TRUE)) {
      tryCatch({
        w_in <- width_mm / 25.4
        h_in <- (height_mm %||% (width_mm * 0.75)) / 25.4
        d <- officer::read_pptx() |>
          officer::add_slide(layout = "Blank", master = "Office Theme")
        d <- rvg::ph_with_vg(d, ggobj = plot, width = w_in, height = h_in)
        pptx_path <- v3_file("figures", paste0(name, "_v3"), "pptx")
        print(d, target = pptx_path)
        message(sprintf("  [pptx] %s", pptx_path))
      }, error = function(e) {
        warning(sprintf("PPTX export failed for %s: %s", name, e$message),
                call. = FALSE)
      })
    }
  }
}

# ── A. 多维多样性时间切片 ─────────────────────────────────────────────

metrics_path <- v3_file("results", paste0("table_community_metrics_with_cri_", RUN_LABEL))
if (file.exists(metrics_path)) {
  m_long <- read_csv(metrics_path, show_col_types = FALSE)
  show_metrics <- c("corrected_richness", "shannon", "pd_prob",
                    "trait_volume", "rao_q")
  m_show <- m_long |>
    filter(metric %in% show_metrics) |>
    mutate(
      metric_label = recode(metric,
        corrected_richness = "Richness",
        shannon = "Shannon",
        pd_prob = "Faith's PD",
        trait_volume = "Trait volume",
        rao_q = "Rao's Q"),
      block_label = factor(block_label, levels = primary_blocks$block_label)
    ) |>
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
    theme(strip.text = element_text(size = 5.5))
  add_plot("fig_v3_multidiversity_timeslices", ts_plot,
           is_map = TRUE, height_mm = 180)
}

# ── B. 趋势 4-panel ──────────────────────────────────────────────────

trend_path <- v3_file("results", paste0("table_grid_trends_with_cri_", RUN_LABEL))
if (file.exists(trend_path)) {
  trd <- read_csv(trend_path, show_col_types = FALSE)
  trd_show <- trd |>
    filter(metric %in% c("corrected_richness", "shannon",
                          "trait_volume", "pd_prob")) |>
    mutate(metric_label = recode(metric,
      corrected_richness = "Richness",
      shannon = "Shannon",
      trait_volume = "Trait volume",
      pd_prob = "Faith's PD")) |>
    group_by(metric_label) |>
    mutate(z = as.numeric(scale(trend_mean)),
           z = pmin(pmax(z, -2.5), 2.5)) |>
    ungroup()
  trd_sf <- grid_sf |> inner_join(trd_show, by = "grid_cell")

  trd_plot <- ggplot(trd_sf) +
    geom_sf(aes(fill = z),
            colour = scales::alpha("white", 0.08), linewidth = 0.04) +
    geom_sf(data = trd_sf |> filter(sig95), aes(geometry = geometry),
            inherit.aes = FALSE, fill = NA,
            colour = "#222222", linewidth = 0.10) +
    geom_sf(data = china_layers, fill = NA, colour = "grey50",
            linewidth = NATURE_LINE_MAP) +
    scale_fill_nature_c(name = "Trend z") +
    facet_wrap(~ metric_label, ncol = 2) +
    coord_sf() +
    theme_nature_map()
  add_plot("fig_v3_community_trends", trd_plot, is_map = TRUE, height_mm = 140)
}

# ── C. varpart vs RF 对比 ─────────────────────────────────────────────

rf_group_path <- v3_file("results", "table_rf_group_importance_summary")
varpart_path  <- v3_file("results",
                          paste0("table_varpart_richness_trend_", RUN_LABEL))
if (file.exists(rf_group_path) && file.exists(varpart_path)) {
  rf_grp <- read_csv(rf_group_path, show_col_types = FALSE)
  vp_df  <- read_csv(varpart_path, show_col_types = FALSE) |>
    filter(grepl("pure", component, ignore.case = TRUE)) |>
    mutate(component = recode(component,
      `Climate (pure)`     = "Climate",
      `Topo+Habitat (pure)` = "Topo+Habitat",
      `Human (pure)`       = "Human",
      `Space (pure)`       = "Space"),
      adj_R2 = as.numeric(adj_R2)) |>
    filter(!is.na(adj_R2))

  compare_plot <- plot_rf_vs_varpart(vp_df, rf_grp)
  add_plot("fig_v3_varpart_vs_rf", compare_plot,
           width_mm = NATURE_WIDTH_L, height_mm = 60)
}

# ── D. RF 变量重要性 ──────────────────────────────────────────────────

rf_var_path <- v3_file("results", "table_rf_importance_summary")
if (file.exists(rf_var_path)) {
  rf_var <- read_csv(rf_var_path, show_col_types = FALSE)
  rf_plot <- plot_rf_variable_importance(rf_var, top_n = 15)
  add_plot("fig_v3_rf_importance", rf_plot,
           width_mm = NATURE_WIDTH_S, height_mm = 80)
}

# ── E. 性状 CWM 空间（含新性状）───────────────────────────────────────

cwm_path <- v3_file("results", paste0("table_cwm_spatial_pattern_", RUN_LABEL))
if (file.exists(cwm_path)) {
  cwm <- read_csv(cwm_path, show_col_types = FALSE)
  cwm_long <- cwm |>
    pivot_longer(matches("^(cwm_|z_)"), names_to = "trait", values_to = "value") |>
    filter(trait %in% c("cwm_diet_specialization", "cwm_habitat_breadth"))
  cwm_sf <- grid_sf |> inner_join(cwm_long, by = "grid_cell")
  if (nrow(cwm_sf) > 0) {
    cwm_plot <- ggplot(cwm_sf) +
      geom_sf(aes(fill = value),
              colour = scales::alpha("white", 0.08), linewidth = 0.04) +
      geom_sf(data = china_layers, fill = NA, colour = "grey50",
              linewidth = NATURE_LINE_MAP) +
      scale_fill_nature_c(name = "CWM", palette = "accent") +
      facet_wrap(~ trait, ncol = 2) +
      coord_sf() +
      theme_nature_map()
    add_plot("fig_v3_cwm_traits", cwm_plot, is_map = TRUE,
             width_mm = NATURE_WIDTH_M, height_mm = 70)
  }
}

# ── F. Master PPTX deck ───────────────────────────────────────────────

if (requireNamespace("officer", quietly = TRUE) &&
    (length(deck_maps) + length(deck_stats)) > 0) {
  log_time("10", "Assembling master PPTX deck...")

  # 地图 deck（PNG 嵌入）
  if (length(deck_maps) > 0) {
    maps_pptx <- officer::read_pptx()
    for (nm in names(deck_maps)) {
      png_path <- v3_file("figures", paste0(nm, "_v3"), "png")
      if (file.exists(png_path)) {
        maps_pptx <- maps_pptx |>
          officer::add_slide(layout = "Blank", master = "Office Theme") |>
          officer::ph_with(
            external_img(png_path, width = 13.33, height = 7.5),
            location = officer::ph_location(width = 13.33, height = 7.5)
          )
      }
    }
    maps_path <- v3_file("figures",
                          paste0("v3_maps_deck_", RUN_LABEL), "pptx")
    print(maps_pptx, target = maps_path)
    message("  [pptx] maps deck → ", maps_path)
  }

  # 统计图 deck（vector dml）
  if (length(deck_stats) > 0 && requireNamespace("rvg", quietly = TRUE)) {
    stats_pptx <- officer::read_pptx()
    for (nm in names(deck_stats)) {
      tryCatch({
        stats_pptx <- stats_pptx |>
          officer::add_slide(layout = "Blank", master = "Office Theme")
        stats_pptx <- rvg::ph_with_vg(stats_pptx, ggobj = deck_stats[[nm]],
                                       width = 12, height = 7)
      }, error = function(e) {
        warning(sprintf("Vector export failed for %s: %s", nm, e$message),
                call. = FALSE)
      })
    }
    stats_path <- v3_file("figures",
                           paste0("v3_stats_deck_editable_", RUN_LABEL), "pptx")
    print(stats_pptx, target = stats_path)
    message("  [pptx] stats deck → ", stats_path)
  }
}

log_time("10", "Done.")
