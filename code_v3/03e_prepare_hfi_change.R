#!/usr/bin/env Rscript
## 03e_prepare_hfi_change.R  —  v3 人类足迹指数变化数据准备
##
## 从本地 HFI 栅格 (data/external/hfi/) 提取各 100km 网格的 HFI 值，
## 计算 delta_hfi (2024 - 2000) 及各 period 的 HFI 值
##
## 数据源: Mu et al. 2022, Scientific Data
##   - 年度全球 HFI, ~1km 分辨率, Mollweide 投影
##   - 本地已有: hfp2000.tif, hfp2005.tif, hfp2010.tif, hfp2015.tif, hfp2020.tif, hfp2024.tif
##
## 注意: HFI 栅格使用 Mollweide 投影，需将网格投影到栅格 CRS 后提取
##   之前的 grid_env 中 HFI 列覆盖率仅 38.2% (476/1247)，需从栅格重新提取
##
## 输出:
##   - data/derived_v3/hfi_change_v3.rds
##   - data/derived_v3/hfi_by_year_v3.rds
##   - results_v3/table_hfi_change_v3.csv

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

log_time("03e", "Starting HFI change data preparation")

# ── 辅助函数：从 HFI 栅格提取各网格均值 ───────────────────────────────
# HFI 栅格为 Mollweide 投影，需将网格投影到栅格 CRS 后提取
extract_hfi_from_raster <- function(tif_path, grid_sf) {
  message(sprintf("    Reading %s (%.0f MB)...",
                  basename(tif_path), file.size(tif_path) / 1e6))
  r <- rast(tif_path)

  # 将网格投影到栅格 CRS (Mollweide)
  grid_proj <- st_transform(grid_sf, crs = crs(r))
  grid_vect_proj <- vect(grid_proj)

  # 裁剪栅格到网格范围（大幅减少内存占用）
  r_crop <- crop(r, ext(grid_vect_proj), snap = "out")

  # 提取各网格均值
  vals <- terra::extract(r_crop, grid_vect_proj, fun = mean, na.rm = TRUE, ID = FALSE)

  col_name <- sprintf("hfi_%s", gsub("\\D", "", tools::file_path_sans_ext(basename(tif_path))))
  tibble(grid_cell = grid_sf$grid_cell, !!col_name := vals[, 1])
}

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

terra::terraOptions(progress = 0)

message(sprintf("[03e] %d grid cells to extract", nrow(grid_sf)))

# ── 2. HFI 栅格路径 ─────────────────────────────────────────────────
hfi_dir <- file.path(PROJECT_ROOT, "data", "external", "hfi")
hfi_years <- c(2000, 2005, 2010, 2015, 2020, 2024)

hfi_paths <- setNames(
  sapply(hfi_years, function(yr) {
    p <- file.path(hfi_dir, sprintf("hfp%d.tif", yr))
    if (file.exists(p)) p else NA_character_
  }),
  hfi_years
)

available_years <- hfi_years[!is.na(hfi_paths)]
# 检查 tif 文件是否有效（非 0 字节）
hfi_paths_valid <- sapply(hfi_paths, function(p) {
  if (is.na(p)) return(NA_character_)
  if (file.size(p) < 1000) return(NA_character_)  # 空文件
  return(p)
})
available_years <- hfi_years[!is.na(hfi_paths_valid)]
message(sprintf("[03e] Available HFI tif years (valid): %s", paste(available_years, collapse = ", ")))

if (length(available_years) < 2) {
  # Fallback: 使用已处理好的 hfi_change_v3.rds
  hfi_existing <- v3_file("derived", "hfi_change_v3", "rds")
  if (file.exists(hfi_existing)) {
    message(sprintf("[03e] No valid HFI tifs. Loading existing: %s", hfi_existing))
    hfi_change <- readRDS(hfi_existing)
    # 过滤到当前站点
    hfi_change <- hfi_change |> filter(grid_cell %in% sites)
    message(sprintf("[03e] Loaded %d grid cells from existing hfi_change_v3.rds", nrow(hfi_change)))
    write_csv(hfi_change, v3_file("results", "table_hfi_change_v3"))
    log_time("03e", sprintf("Completed (from existing rds, %d cells)", nrow(hfi_change)))
    quit(save = "no", status = 0)
  }
  stop("[03e] Need at least 2 HFI raster years and no existing rds fallback! Available: ",
       paste(available_years, collapse = ", "))
}

# ── 3. 从栅格提取各年份 HFI ──────────────────────────────────────────
# 直接从栅格提取（grid_env 中的 HFI 覆盖率仅 38%，不可靠）
# 策略：先提取较小的 2020/2024 栅格（~570MB），再提取 2000-2015（~2.4GB）

message("[03e] Extracting HFI from rasters (projecting grid to Mollweide CRS)...")

# 按文件大小排序：小文件先提取（快速验证流程）
yr_order <- available_years[order(sapply(hfi_paths[as.character(available_years)], function(p) {
  if (is.na(p)) Inf else file.size(p)
}))]

hfi_list <- list()
for (yr in yr_order) {
  tif_path <- hfi_paths[as.character(yr)]
  if (is.na(tif_path) || !file.exists(tif_path)) next

  hfi_yr <- extract_hfi_from_raster(tif_path, grid_sf)
  n_valid <- sum(!is.na(hfi_yr[[2]]))
  message(sprintf("[03e]   Year %d: %d/%d grids with valid data (%.1f%%)",
                  yr, n_valid, nrow(grid_sf), 100 * n_valid / nrow(grid_sf)))
  hfi_list[[as.character(yr)]] <- hfi_yr

  # 释放内存
  gc()
}

# 合并所有年份
hfi_all <- Reduce(function(df1, df2) left_join(df1, df2, by = "grid_cell"), hfi_list)

available_hfi_cols <- grep("^hfi_\\d{4}$", names(hfi_all), value = TRUE)
available_hfi_years <- as.integer(gsub("hfi_", "", available_hfi_cols))

message(sprintf("[03e] Extracted HFI years: %s", paste(available_hfi_years, collapse = ", ")))

# 检查覆盖率
for (col in available_hfi_cols) {
  n_valid <- sum(!is.na(hfi_all[[col]]))
  message(sprintf("[03e]   %s: %d/%d valid (%.1f%%)",
                  col, n_valid, nrow(hfi_all), 100 * n_valid / nrow(hfi_all)))
}

# ── 4. 用栅格提取值更新 grid_environment ─────────────────────────
grid_env_stem <- paste0("grid_environment", GRID_TAG, "_v3")
grid_env_path <- v3_file("derived", grid_env_stem, "rds")
if (!file.exists(grid_env_path)) grid_env_path <- v3_file("derived", "grid_environment_v3", "rds")
if (file.exists(grid_env_path)) {
  grid_env <- readRDS(grid_env_path)
  old_hfi_cols <- grep("^hfi_", names(grid_env), value = TRUE)
  message(sprintf("[03e] Updating grid_environment_v3 HFI columns (was %d cols, %d valid)",
                  length(old_hfi_cols),
                  sum(!is.na(grid_env$hfi_2000))))

  # 移除旧 HFI 列，替换为新提取值
  grid_env <- grid_env |> select(-any_of(old_hfi_cols))
  grid_env <- grid_env |> left_join(hfi_all, by = "grid_cell")

  # 重新计算 hfi_mean 和 hfi_sd
  grid_env <- grid_env |>
    rowwise() |>
    mutate(
      hfi_mean = mean(c_across(all_of(available_hfi_cols)), na.rm = TRUE),
      hfi_sd   = sd(c_across(all_of(available_hfi_cols)), na.rm = TRUE)
    ) |>
    ungroup()

  saveRDS(grid_env, grid_env_path)
  message(sprintf("[03e] grid_environment_v3 updated: hfi_2000 valid = %d",
                  sum(!is.na(grid_env$hfi_2000))))
}

# ── 5. 计算 HFI 变化量 ──────────────────────────────────────────────
early_yr <- min(available_hfi_years)
late_yr  <- max(available_hfi_years)

early_col <- sprintf("hfi_%d", early_yr)
late_col  <- sprintf("hfi_%d", late_yr)

message(sprintf("[03e] Computing HFI change: %d -> %d", early_yr, late_yr))

hfi_change <- hfi_all |>
  mutate(
    delta_hfi = .data[[late_col]] - .data[[early_col]],
    early_year = early_yr,
    late_year  = late_yr,
    data_source = "HFI_Mu_2022"
  ) |>
  left_join(
    grid_sf |>
      st_drop_geometry() |>
      as_tibble() |>
      select(grid_cell, centroid_lon, centroid_lat),
    by = "grid_cell"
  ) |>
  relocate(grid_cell, centroid_lon, centroid_lat)

# ── 6. 计算各 period 的 HFI 均值 ────────────────────────────────────
# P1 (2000-2004): hfi_2000
# P2 (2005-2009): hfi_2005
# P3 (2010-2014): hfi_2010
# P4 (2015-2019): hfi_2015
# P5 (2020-2024): (hfi_2020 + hfi_2024) / 2

period_hfi <- hfi_all |>
  select(grid_cell, all_of(available_hfi_cols)) |>
  pivot_longer(-grid_cell, names_to = "year_str", values_to = "hfi") |>
  mutate(year = as.integer(gsub("hfi_", "", year_str))) |>
  select(-year_str)

period_map <- tibble(
  year = available_hfi_years,
  period = case_when(
    year <= 2004 ~ "P1", year <= 2009 ~ "P2", year <= 2014 ~ "P3",
    year <= 2019 ~ "P4", TRUE ~ "P5"
  )
)

period_hfi <- period_hfi |>
  left_join(period_map, by = "year") |>
  group_by(grid_cell, period) |>
  summarise(hfi_period = mean(hfi, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = period, values_from = hfi_period,
              names_prefix = "hfi_")

hfi_change <- hfi_change |>
  left_join(period_hfi, by = "grid_cell")

# ── 7. 保存 ──────────────────────────────────────────────────────────
saveRDS(hfi_change, v3_file("derived", "hfi_change_v3", "rds"))
saveRDS(hfi_all, v3_file("derived", "hfi_by_year_v3", "rds"))
write_csv(hfi_change, v3_file("results", "table_hfi_change_v3"))

message(sprintf("[03e] %d grid cells with HFI change data", nrow(hfi_change)))
message(sprintf("[03e] Period: %d -> %d", early_yr, late_yr))
message(sprintf("[03e] delta_hfi range: [%.2f, %.2f] (mean %.2f)",
                min(hfi_change$delta_hfi, na.rm = TRUE),
                max(hfi_change$delta_hfi, na.rm = TRUE),
                mean(hfi_change$delta_hfi, na.rm = TRUE)))
message(sprintf("[03e] N valid delta_hfi: %d / %d",
                sum(!is.na(hfi_change$delta_hfi)), nrow(hfi_change)))

log_time("03e", "Completed")
