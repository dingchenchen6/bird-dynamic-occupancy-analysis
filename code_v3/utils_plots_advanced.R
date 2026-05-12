#!/usr/bin/env Rscript
## utils_plots_advanced.R  —  v3 高级绘图工具
##
## Forest plot、raincloud plot、Nature 风格扩展
## 修复 v2 中 %||% 重复定义问题（统一使用 utils_paths 中的版本）

# ── 加载依赖 ──────────────────────────────────────────────────────────
{
  .root <- Sys.getenv("BIRD_PROJECT_ROOT",
    if (dir.exists(file.path("~", "Documents", "New project",
                             "bird_dynamic_occupancy_analysis")))
      file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis")
    else getwd())
  source(file.path(.root, "code_v3", "utils_paths.R"))
  rm(.root)
}

# ── Forest Plot（Nature 风格）─────────────────────────────────────────
#' 绘制回归系数 forest plot，带 95% CRI
#'
#' @param coef_df tibble: term, estimate, q025, q975, group (可选)
#' @param title 标题
#' @param vline 竖线位置（默认 0）
#' @param color_by_group 是否按组着色
forest_plot_nature <- function(coef_df, title = NULL, vline = 0,
                                color_by_group = FALSE) {
  p <- ggplot(coef_df, aes(y = term, x = estimate)) +
    geom_vline(xintercept = vline, linetype = "dashed",
               colour = "grey50", linewidth = 0.2) +
    geom_linerange(aes(xmin = q025, xmax = q975),
                   linewidth = 0.4, linewidth_head = 0.8) +
    geom_point(size = 1.2, shape = 18)

  if (color_by_group && "group" %in% names(coef_df)) {
    p <- p + aes(colour = group) +
      scale_colour_manual(values = pal_nature_discrete(length(unique(coef_df$group)))) +
      guides(colour = guide_legend(override.aes = list(size = 1.5)))
  }

  p +
    labs(title = title, x = "Effect size (95% CRI)", y = NULL) +
    theme_nature_pub() +
    theme(axis.text.y = element_text(face = "italic"))
}

# ── Raincloud Plot ────────────────────────────────────────────────────
raincloud_nature <- function(df, x_var, y_var, fill_var = NULL,
                              title = NULL) {
  x_sym <- rlang::ensym(x_var)
  y_sym <- rlang::ensym(y_var)

  p <- ggplot(df, aes(x = {{ x_var }}, y = {{ y_var }}))

  if (!is.null(fill_var)) {
    p <- p + aes(fill = {{ fill_var }})
  }

  # Half-violin
  p <- p +
    gghalves::geom_half_violin(side = "l", adjust = 0.6,
                                width = 0.5, alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.5) +
    geom_jitter(size = 0.3, alpha = 0.2, width = 0.05) +
    labs(title = title) +
    theme_nature_pub()

  p
}

# ── 多面板趋势轨迹图 ──────────────────────────────────────────────────
trajectory_plot_nature <- function(trend_df, metric_col = "metric",
                                    value_col = "mean",
                                    ci_lo = "q025", ci_hi = "q975",
                                    facet_var = NULL) {
  p <- ggplot(trend_df,
              aes(x = period, y = .data[[value_col]])) +
    geom_ribbon(aes(ymin = .data[[ci_lo]], ymax = .data[[ci_hi]]),
                alpha = 0.15, fill = NATURE_ACCENT) +
    geom_line(colour = NATURE_ACCENT, linewidth = 0.5) +
    geom_point(colour = NATURE_ACCENT, size = 1.2) +
    labs(x = "Period", y = "Value") +
    theme_nature_pub()

  if (!is.null(facet_var)) {
    p <- p + facet_wrap(vars(.data[[facet_var]]), scales = "free_y")
  }

  p
}

# ── 配色辅助 ──────────────────────────────────────────────────────────
# FIX #12: pal_nature_discrete 统一定义在 utils_mapping.R 中
# 此处不再重复定义，避免 source 顺序不同导致版本分叉
# 若 utils_mapping 未加载则提供 fallback
if (!exists("pal_nature_discrete", mode = "function")) {
  pal_nature_discrete <- function(n = 8) {
    cols <- c("#0E5A78", "#B8860B", "#6B8E23", "#CD5C5C",
              "#708090", "#DAA520", "#8B4513", "#4682B4")
    if (n <= length(cols)) cols[seq_len(n)] else rep_len(cols, n)
  }
}

message("[utils_plots_advanced] loaded")
