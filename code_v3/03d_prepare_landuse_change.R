#!/usr/bin/env Rscript
## 03d_prepare_landuse_change.R  —  v3 土地利用变化数据准备
##
## 使用 CLCD (China Land Cover Dataset, Yang & Huang 2026)
## 提取各 100km 网格内年度土地覆盖类型面积比例，
## 计算 P1→P5 (2000→2020) 的土地利用变化量
##
## 数据源: https://doi.org/10.5281/zenodo.18180184
##   - CLCD: 30m 分辨率, 年度 1985-2025, 9 类地覆盖
##   - 投影: Albers Equal Area
##
## 输出:
##   - data/derived_v3/landuse_change_v3.rds
##   - results_v3/table_landuse_change_v3.csv

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(tidyr); library(tibble)
  library(readr)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
P <- ensure_v3_dirs()

log_time("03d", "Starting land use change data preparation")

# ── 1. 加载网格 ──────────────────────────────────────────────────────
grid_rds_path <- v3_file("derived", paste0("china_grid_", GRID_SIZE_KM, "km_v3"), "rds")
grid_sf <- if (file.exists(grid_rds_path)) {
  readRDS(grid_rds_path)
} else {
  safe_read(file.path(DIRS$v2_derived, "china_grid_100km_v2.rds"))
}
survey_stem <- paste0("survey_history", GRID_TAG, "_v3")
survey_path <- v3_file("derived", survey_stem, "rds")
if (!file.exists(survey_path)) survey_path <- v3_file("derived", "survey_history_v3", "rds")
sites <- readRDS(survey_path)$sites
grid_sf <- grid_sf |> filter(grid_cell %in% sites)

# CLCD 使用 Albers Equal Area 投影，需要将网格转换
clcd_crs <- "+proj=aea +lat_1=25 +lat_2=47 +lat_0=0 +lon_0=105 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
grid_albers <- st_transform(grid_sf, crs = clcd_crs)
grid_vect_albers <- vect(grid_albers)

# ── 2. CLCD 数据路径 ──────────────────────────────────────────────────
clcd_dir <- Sys.getenv("V3_CLCD_DIR",
  file.path(PROJECT_ROOT, "data", "external", "clcd"))
if (!dir.exists(clcd_dir)) dir.create(clcd_dir, recursive = TRUE, showWarnings = FALSE)

# CLCD 文件命名: CLCD_v01_YYYY_albert.tif
# 需要的关键年份: 2000 (P1), 2005 (P2), 2010 (P3), 2015 (P4), 2020 (P5)
# CLCD v1.0.5 覆盖到 2025 (DOI: 10.5281/zenodo.18180184)
# 2024 为可选：有则用，无则 fallback 到 2020 作为最新时间点
clcd_years <- c(2000, 2005, 2010, 2015, 2020, 2024)

# CLCD 分类编码:
# 1 = Cropland, 2 = Forest, 3 = Shrub, 4 = Grassland,
# 5 = Water, 6 = Snow/Ice, 7 = Barren, 8 = Impervious, 9 = Wetland
CLCD_CLASSES <- c(
  "cropland" = 1, "forest" = 2, "shrub" = 3, "grassland" = 4,
  "water" = 5, "snow_ice" = 6, "barren" = 7, "impervious" = 8, "wetland" = 9
)

# 生态意义分组
NATURAL_CLASSES <- c("forest", "shrub", "grassland", "water", "wetland")

# ── 3. 检查/下载 CLCD 数据 ──────────────────────────────────────────
clcd_tif_paths <- sapply(clcd_years, function(yr) {
  # 尝试多种命名模式
  candidates <- c(
    file.path(clcd_dir, sprintf("CLCD_v01_%d_albert.tif", yr)),
    file.path(clcd_dir, sprintf("CLCD_v01_%d_albert_province", yr)),
    file.path(clcd_dir, sprintf("clcd_%d.tif", yr))
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) > 0) return(found[1])
  return(NA_character_)
})

missing_years <- clcd_years[is.na(clcd_tif_paths)]

if (length(missing_years) > 0) {
  message(sprintf("[03d] CLCD tif files not found for years: %s",
                  paste(missing_years, collapse = ", ")))
  message(sprintf("[03d] Expected location: %s", clcd_dir))
  message("[03d] Please download CLCD from https://doi.org/10.5281/zenodo.18180184")
  message("[03d] The full dataset is ~60GB. You only need these files:")
  for (yr in missing_years) {
    message(sprintf("  - CLCD_v01_%d_albert.tif", yr))
  }
  message("[03d] Each file is ~2GB.")

  # 备选方案：如果已下载部分年份，用已有数据
  available_years <- clcd_years[!is.na(clcd_tif_paths)]
  if (length(available_years) == 0) {
    message("[03d] No CLCD data available. Using fallback from existing landcover.")

    # 备选：从 grid_env 的现有 landcover 数据构建（仅单时间点，无法计算变化）
    grid_env_stem <- paste0("grid_environment", GRID_TAG, "_v3")
    grid_env_path <- v3_file("derived", grid_env_stem, "rds")
    if (!file.exists(grid_env_path)) grid_env_path <- v3_file("derived", "grid_environment_v3", "rds")
    grid_env <- readRDS(grid_env_path)

    # 使用 HFI 变化作为人类压力变化的代理
    # landcover 变化用占位（均为 0 变化，即假设不变）
    landuse_change <- grid_env |>
      select(grid_cell, centroid_lon, centroid_lat) |>
      mutate(
        delta_forest     = 0,
        delta_shrub      = 0,
        delta_grassland  = 0,
        delta_wetland    = 0,
        delta_water      = 0,
        delta_impervious = 0,
        delta_cropland   = 0,
        delta_natural    = 0,
        data_source      = "placeholder"
      )

    saveRDS(landuse_change, v3_file("derived", "landuse_change_v3", "rds"))
    write_csv(landuse_change, v3_file("results", "table_landuse_change_v3"))
    warning("[03d] Using PLACEHOLDER landuse deltas. Download CLCD for real values!")
    log_time("03d", "Completed (PLACEHOLDER — need CLCD)")
    quit(save = "no", status = 0)
  }

  # 有部分年份可用
  clcd_years <- available_years
  clcd_tif_paths <- clcd_tif_paths[!is.na(clcd_tif_paths)]
}

# ── 4. 提取各年份各网格的土地覆盖面积比例 ────────────────────────────
# 这是核心步骤：对每个 CLCD 栅格，统计每个 100km 网格内各类的面积比例

extract_landcover_proportions <- function(tif_path, grid_vect, class_codes) {
  message(sprintf("[03d] Reading %s", basename(tif_path)))
  r <- rast(tif_path)

  # 确保投影一致
  if (crs(r) != crs(grid_vect)) {
    grid_vect <- project(grid_vect, crs(r))
  }

  n_cells <- length(grid_vect)
  n_classes <- length(class_codes)
  result <- matrix(NA_real_, nrow = n_cells, ncol = n_classes)
  colnames(result) <- names(class_codes)

  message(sprintf("[03d] Extracting proportions for %d grid cells, %d classes...",
                  n_cells, n_classes))

  # 方法1: exactextractr（推荐，更精确）
  if (requireNamespace("exactextractr", quietly = TRUE)) {
    grid_sf_albers <- st_as_sf(grid_vect)
    for (i in seq_along(class_codes)) {
      class_name <- names(class_codes)[i]
      class_val <- class_codes[i]

      # 创建二值栅格：1 = 该类，0/NA = 其他
      r_binary <- ifel(r == class_val, 1, 0)
      props <- exactextractr::exact_extract(r_binary, grid_sf_albers, "mean")
      result[, i] <- props
    }
  } else {
    # 方法2: terra::extract（较快但精度稍低）
    for (i in seq_along(class_codes)) {
      class_name <- names(class_codes)[i]
      class_val <- class_codes[i]

      r_binary <- ifel(r == class_val, 1, 0)
      # 对 30m 栅格用 100km 网格聚合 — 先降低分辨率
      r_agg <- aggregate(r_binary, fact = 10, fun = "mean", na.rm = TRUE)
      ext_vals <- terra::extract(r_agg, grid_vect, fun = "mean", na.rm = TRUE)
      result[, i] <- ext_vals[, 2]
    }
  }

  as_tibble(result) |>
    mutate(grid_cell = grid_sf$grid_cell) |>
    relocate(grid_cell)
}

# 提取各年份
landcover_by_year <- list()

for (yr in clcd_years) {
  tif_path <- clcd_tif_paths[as.character(yr)]
  if (is.na(tif_path) || !file.exists(tif_path)) next

  lc <- extract_landcover_proportions(tif_path, grid_vect_albers, CLCD_CLASSES)

  # 重命名列加年份后缀
  lc_renamed <- lc |>
    rename_with(~ paste0(.x, "_", yr), -grid_cell)

  landcover_by_year[[as.character(yr)]] <- lc_renamed
  message(sprintf("[03d] Year %d: %d grids extracted", yr, nrow(lc_renamed)))
}

# 合并所有年份
if (length(landcover_by_year) == 0) {
  stop("[03d] No land cover data extracted!")
}

lc_all <- landcover_by_year[[1]]
for (i in seq_along(landcover_by_year)[-1]) {
  lc_all <- lc_all |>
    left_join(landcover_by_year[[i]], by = "grid_cell")
}

# ── 5. 计算土地利用变化量 ──────────────────────────────────────────
# 使用最早和最晚可用年份
early_year <- min(clcd_years)
late_year <- max(clcd_years)

message(sprintf("[03d] Computing land use change: %d -> %d", early_year, late_year))

# 需要计算 delta 的类型
lc_types <- c("forest", "shrub", "grassland", "wetland", "water", "impervious", "cropland")

landuse_change <- tibble(grid_cell = lc_all$grid_cell)

for (lt in lc_types) {
  early_col <- paste0(lt, "_", early_year)
  late_col  <- paste0(lt, "_", late_year)

  if (early_col %in% names(lc_all) && late_col %in% names(lc_all)) {
    landuse_change[[paste0("delta_", lt)]] <- lc_all[[late_col]] - lc_all[[early_col]]
  } else {
    landuse_change[[paste0("delta_", lt)]] <- NA_real_
  }
}

# 自然地类总比例变化
natural_early <- rowSums(select(lc_all, intersect(
  paste0(NATURAL_CLASSES, "_", early_year), names(lc_all))), na.rm = TRUE)
natural_late <- rowSums(select(lc_all, intersect(
  paste0(NATURAL_CLASSES, "_", late_year), names(lc_all))), na.rm = TRUE)
landuse_change$delta_natural <- natural_late - natural_early

# 加上元数据
landuse_change <- landuse_change |>
  left_join(
    grid_sf |>
      st_drop_geometry() |>
      as_tibble() |>
      select(grid_cell, centroid_lon, centroid_lat),
    by = "grid_cell"
  ) |>
  relocate(grid_cell, centroid_lon, centroid_lat) |>
  mutate(
    early_year = early_year,
    late_year  = late_year,
    data_source = "CLCD"
  )

# ── 6. 中间 period 的 landcover 也保存（用于 period-level 分析）───
# 将完整的 lc_all 也保存（含各年份面积比例）
landuse_full <- lc_all |>
  left_join(
    grid_sf |>
      st_drop_geometry() |>
      as_tibble() |>
      select(grid_cell, centroid_lon, centroid_lat),
    by = "grid_cell"
  ) |>
  relocate(grid_cell, centroid_lon, centroid_lat)

# ── 7. 保存 ──────────────────────────────────────────────────────────
saveRDS(landuse_change, v3_file("derived", "landuse_change_v3", "rds"))
saveRDS(landuse_full, v3_file("derived", "landuse_by_year_v3", "rds"))
write_csv(landuse_change, v3_file("results", "table_landuse_change_v3"))

message(sprintf("[03d] %d grid cells with landuse change data", nrow(landuse_change)))
message(sprintf("[03d] Period: %d -> %d", early_year, late_year))
for (lt in lc_types) {
  delta_col <- paste0("delta_", lt)
  if (delta_col %in% names(landuse_change)) {
    message(sprintf("[03d]   %s: [%.4f, %.4f] (mean %.4f)",
                    lt,
                    min(landuse_change[[delta_col]], na.rm = TRUE),
                    max(landuse_change[[delta_col]], na.rm = TRUE),
                    mean(landuse_change[[delta_col]], na.rm = TRUE)))
  }
}

log_time("03d", "Completed")
