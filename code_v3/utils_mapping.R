#!/usr/bin/env Rscript
## utils_mapping.R  —  v3 制图工具模块
##
## Nature/Science 风格地图和出版级图表主题
## 修复 v2 中 %||% 在 L270 定义却在 L264 使用的 source-order bug

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

# ── Nature 出版级主题（非地图）────────────────────────────────────────
theme_nature_pub <- function(base_size = 7.5, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # 面板
      panel.border      = element_rect(colour = "black", linewidth = 0.25),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.background  = element_rect(fill = "white", colour = NA),

      # 坐标轴
      axis.line         = element_line(colour = "black", linewidth = 0.25),
      axis.ticks        = element_line(colour = "black", linewidth = 0.15),
      axis.text         = element_text(colour = "black", size = base_size),
      axis.title        = element_text(colour = "black", size = base_size,
                                        face = "plain"),

      # 图例
      legend.background  = element_rect(fill = "white", colour = NA),
      legend.key         = element_rect(fill = "white", colour = NA),
      legend.key.size    = unit(3, "pt"),
      legend.text        = element_text(size = base_size - 1),
      legend.title       = element_text(size = base_size, face = "plain"),
      legend.position    = "right",

      # strip（facet）
      strip.background  = element_rect(fill = "grey92", colour = "black",
                                        linewidth = 0.15),
      strip.text        = element_text(size = base_size, face = "plain",
                                        margin = margin(2, 2, 2, 2)),

      # 标题
      plot.title        = element_text(size = base_size + 1, face = "bold",
                                        hjust = 0),
      plot.subtitle     = element_text(size = base_size, hjust = 0,
                                        colour = "grey30"),
      plot.margin       = margin(4, 4, 4, 4, unit = "pt")
    )
}

# ── Nature 地图主题 ──────────────────────────────────────────────────
theme_nature_map <- function(base_size = 7.5, base_family = "Arial") {
  theme_nature_pub(base_size, base_family) %+replace%
    theme(
      axis.line       = element_blank(),
      axis.ticks      = element_blank(),
      axis.text       = element_text(size = base_size - 1.5, colour = "grey40"),
      axis.title      = element_blank(),
      panel.border    = element_rect(colour = "grey60", linewidth = 0.15),
      panel.background = element_rect(fill = "white", colour = NA),
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.key.width = unit(6, "mm"),
      legend.key.height = unit(2, "mm")
    )
}

# ── 保存 Nature 级图表 ──────────────────────────────────────────────
#' @param gg ggplot 对象
#' @param stem 文件名（无扩展名）
#' @param width_mm 宽度（mm），默认 89（单栏）
#' @param height_mm 高度（mm）
#' @param dpi 分辨率
save_nature <- function(gg, stem, width_mm = 89, height_mm = NULL,
                         dpi = 300, dir_sub = "figures") {
  load_config_if_needed()
  out_dir <- DIRS[[dir_sub]] %||% DIRS$figures
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  width_in  <- width_mm / 25.4
  # 自动高度
  if (is.null(height_mm)) {
    # 从 gg 构建：获取 aspect ratio
    height_in <- width_in * 0.75  # 默认 4:3
  } else {
    height_in <- height_mm / 25.4
  }

  # PNG
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  ggsave(png_path, gg, width = width_in, height = height_in,
         dpi = dpi, bg = "white")
  # PDF
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  ggsave(pdf_path, gg, width = width_in, height = height_in,
         device = cairo_pdf, bg = "white")

  message(sprintf("[save_nature] %s (%.0f×%.0fmm) → png+pdf",
                  stem, width_mm, height_mm %||% (width_mm * 0.75)))
  invisible(c(png = png_path, pdf = pdf_path))
}

# ── 配色方案 ─────────────────────────────────────────────────────────

# 灰度 + 单色点缀
pal_nature_accent <- function(n = 5) {
  c("#0E5A78", "#4A7A96", "#7A9AAE", "#AABAC6", "#D4DEE4")[seq_len(n)]
}

# 散点双向色带（Blue-White-Red）
pal_div_bwr <- function(n = 100) {
  colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(n)
}

# 散点 LaJolla 色带
pal_lajolla <- function(n = 100) {
  if (!requireNamespace("scico", quietly = TRUE)) {
    pal_div_bwr(n)
  } else {
    scico::scico(n, palette = "lajolla")
  }
}

# 离散分类色板
pal_nature_discrete <- function(n = 8) {
  cols <- c("#0E5A78", "#B8860B", "#6B8E23", "#CD5C5C",
            "#708090", "#DAA520", "#8B4513", "#4682B4")
  if (n <= length(cols)) cols[seq_len(n)] else rep_len(cols, n)
}

# ── 中国省界底图 ──────────────────────────────────────────────────────
china_province_basemap <- function() {
  load_config_if_needed()
  shp_path <- file.path(DIRS$data_raw, "china_province.shp")
  if (file.exists(shp_path)) {
    sf::st_read(shp_path, quiet = TRUE)
  } else {
    # 备选：rnaturalearth
    if (requireNamespace("rnaturalearth", quietly = TRUE)) {
      rnaturalearth::ne_countries(country = "china", scale = 50,
                                   returnclass = "sf")
    } else {
      warning("No China boundary shapefile found at ", shp_path)
      NULL
    }
  }
}

# ── 连续色标（Nature 风格）────────────────────────────────────────────
scale_fill_nature_c <- function(name = waiver(), palette = "accent", ...) {
  if (palette == "accent") {
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0, name = name, ...)
  } else if (palette == "lajolla") {
    if (requireNamespace("scico", quietly = TRUE)) {
      scico::scale_fill_scico(palette = "lajolla", name = name, ...)
    } else {
      scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                           midpoint = 0, name = name, ...)
    }
  } else {
    scale_fill_gradient(name = name, ...)
  }
}

message("[utils_mapping] loaded")
