#!/usr/bin/env Rscript
## 03_prepare_environment.R  —  v3 环境变量准备
##
## 提取并标准化环境协变量，统一年份范围
## FIX: 10km 网格不再依赖 v2 的 100km grid_environment
##      改为从 10km 网格质心提取环境变量（与 02 生成的 grid_cell 对齐）

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(sf)
  library(terra)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_spatial.R"))
P <- ensure_v3_dirs()

log_time("03", "Starting environment preparation")

# ── 1. 加载网格和站点信息 ──────────────────────────────────────────────
survey_stem <- paste0("survey_history", GRID_TAG, "_v3")
survey <- safe_read(v3_file("derived", survey_stem, "rds"))
if (is.null(survey)) {
  survey <- safe_read(v3_file("derived", "survey_history_v3", "rds"))
}
if (is.null(survey)) {
  survey <- safe_read(file.path(DIRS$v2_derived, "survey_history_v2.rds"))
}
sites <- survey$sites
n_sites <- length(sites)
message(sprintf("[03] %d sites from survey (%dkm grid)", n_sites, GRID_SIZE_KM))

# 加载与 02 对齐的网格 sf
grid_rds_path <- v3_file("derived", paste0("china_grid_", GRID_SIZE_KM, "km_v3"), "rds")
grid_sf <- safe_read(grid_rds_path)

if (is.null(grid_sf)) {
  # fallback: 加载 v2 100km 网格
  grid_sf <- safe_read(file.path(DIRS$v2_derived, "china_grid_100km_v2.rds"))
}

# ── 2. 加载环境数据 ──────────────────────────────────────────────────
# 策略：优先尝试与当前网格对齐的 v3 环境数据
# 如果是 10km 网格，需要从栅格重新提取（v2 是 100km 的）

grid_env <- NULL

# 尝试 1: v2 grid_environment（仅适用于 100km 网格）
if (GRID_SIZE_KM == 100L) {
  grid_env <- safe_read(file.path(DIRS$v2_derived,
                                    "grid_environment_dynamic_occupancy.rds"))
  if (!is.null(grid_env)) {
    # 过滤到当前 survey 的站点
    grid_env <- grid_env |>
      filter(grid_cell %in% sites) |>
      arrange(match(grid_cell, sites))
    n_matched <- nrow(grid_env)
    message(sprintf("[03] v2 grid_env: %d / %d sites matched", n_matched, n_sites))
    if (n_matched < n_sites * 0.5) {
      warning("[03] v2 grid_env matches <50% of sites. Falling back to raster extraction.")
      grid_env <- NULL
    }
  }
}

# 尝试 2: 同分辨率的已有 v3 grid_environment
grid_env_stem <- paste0("grid_environment", GRID_TAG, "_v3")
if (is.null(grid_env)) {
  grid_env <- safe_read(v3_file("derived", grid_env_stem, "rds"))
  if (!is.null(grid_env)) {
    grid_env <- grid_env |> filter(grid_cell %in% sites)
    message(sprintf("[03] v3 grid_env: %d / %d sites matched", nrow(grid_env), n_sites))
  }
}

# 尝试 3: 从栅格数据重新提取（适用于 10km 或 grid_env 匹配不足的情况）
if (is.null(grid_env) || nrow(grid_env) < n_sites * 0.5) {
  message(sprintf("[03] Extracting environment from raster for %dkm grid (%d sites)",
                  GRID_SIZE_KM, n_sites))
  grid_env <- extract_env_from_raster(grid_sf, sites)
}

if (is.null(grid_env) || nrow(grid_env) == 0) {
  stop("[03] Cannot prepare environment data. Check grid and raster availability.")
}

# 确保关键列存在
env_cols <- c("bio4", "bio7", "bio11", "bio13",
               "elev_mean", "elev_sd", "texture_shannon",
               "habitat_diversity_shannon",
               "hfi_mean", "landcover_built", "landcover_cropland",
               "centroid_lon", "centroid_lat")

missing_cols <- setdiff(env_cols, names(grid_env))
if (length(missing_cols) > 0) {
  warning("[03] Missing environment columns: ", paste(missing_cols, collapse = ", "))
  for (mc in missing_cols) grid_env[[mc]] <- NA_real_
}

# year_scaled（period-level 中心化，在模型调用时动态设置）
grid_env$year_scaled <- 0

# ── 3. 检测协变量 ──────────────────────────────────────────────────────
# 从 survey$visits（宽格式 grid_cell × P1:P5）构建
visits_wide <- survey$visits
period_cols <- grep("^P[0-9]+$", names(visits_wide), value = TRUE)
det_cov <- visits_wide |>
  pivot_longer(all_of(period_cols), names_to = "period", values_to = "n_events") |>
  filter(grid_cell %in% sites) |>
  mutate(log_events = log1p(n_events))

# 尝试从 v2 visit_effort 提取 duration 信息补充
visit_effort <- safe_read(file.path(DIRS$v2_derived, "visit_effort_2000_2024.rds"))
if (!is.null(visit_effort) && "block_id" %in% names(visit_effort)) {
  # FIX: 仅 merge 当前 sites 范围内的 visit_effort
  # 10km grid_cell 与 v2 grid_cell 编码不同，按 grid_cell 数值做交集
  ve_sites <- intersect(visit_effort$grid_cell, sites)
  if (length(ve_sites) > n_sites * 0.3) {
    # 有足够交集时才使用 v2 duration
    duration_by_cell_period <- visit_effort |>
      filter(grid_cell %in% ve_sites) |>
      mutate(period = paste0("P", block_id)) |>
      group_by(grid_cell, period) |>
      summarise(
        log_duration_mean = mean(log10(pmax(mean_duration_min, 1)), na.rm = TRUE),
        has_duration = as.integer(any(!is.na(mean_duration_min) & mean_duration_min > 0)),
        .groups = "drop"
      ) |>
      mutate(log_duration = ifelse(has_duration == 1, log_duration_mean, 0))
    det_cov <- det_cov |>
      left_join(duration_by_cell_period, by = c("grid_cell", "period")) |>
      mutate(
        has_duration = replace_na(has_duration, 0L),
        log_duration = replace_na(log_duration, 0)
      )
  } else {
    message("[03] v2 visit_effort has low overlap with current sites. Using default detection covariates.")
    det_cov$has_duration <- 0L
    det_cov$log_duration <- 0
  }
} else {
  det_cov$has_duration <- 0L
  det_cov$log_duration <- 0
}

# ── 4. 保存 ──────────────────────────────────────────────────────────
grid_env_stem <- paste0("grid_environment", GRID_TAG, "_v3")
det_cov_stem  <- paste0("detection_covariates", GRID_TAG, "_v3")
saveRDS(grid_env, v3_file("derived", grid_env_stem, "rds"))
saveRDS(det_cov, v3_file("derived", det_cov_stem, "rds"))

write_csv(grid_env, v3_file("results", paste0("table_grid_environment", GRID_TAG, "_v3")))
message(sprintf("[03] %d grid cells × %d environment variables",
                nrow(grid_env), length(env_cols)))

log_time("03", sprintf("Completed: %d sites matched", nrow(grid_env)))

# ── 辅助函数：从栅格提取环境变量 ──────────────────────────────────────
extract_env_from_raster <- function(grid_sf, sites) {
  # 从网格 sf 中提取质心坐标，再从栅格提取环境变量
  if (is.null(grid_sf)) {
    warning("[03] grid_sf not available for raster extraction")
    return(NULL)
  }

  # 筛选当前 sites 的网格
  grid_subset <- grid_sf |> filter(grid_cell %in% sites)
  if (nrow(grid_subset) == 0) {
    warning("[03] No grid cells match current sites")
    return(NULL)
  }

  # 计算质心
  centroids <- grid_subset |>
    st_centroid() |>
    mutate(
      centroid_lon = st_coordinates(geometry)[, 1],
      centroid_lat = st_coordinates(geometry)[, 2]
    ) |>
    st_drop_geometry() |>
    select(grid_cell, centroid_lon, centroid_lat)

  # 尝试加载 WorldClim / CHELSA 生物气候栅格
  bio_vars <- c("bio4", "bio7", "bio11", "bio13")
  bio_stack <- NULL
  bio_dirs <- c(
    file.path(DIRS$data_raw, "worldclim"),
    file.path(DIRS$data_raw, "chelsa"),
    file.path(DIRS$data_raw, "external", "worldclim")
  )
  for (d in bio_dirs) {
    if (dir.exists(d)) {
      bio_files <- list.files(d, pattern = "^bio[0-9]+.*\\.tif$", full.names = TRUE)
      if (length(bio_files) >= 4) {
        bio_stack <- tryCatch(rast(bio_files), error = function(e) NULL)
        if (!is.null(bio_stack)) {
          message(sprintf("[03] Loaded bioclimate from %s (%d layers)", d, nlyr(bio_stack)))
          break
        }
      }
    }
  }

  # 尝试加载 DEM
  elev_rast <- NULL
  elev_paths <- c(
    file.path(DIRS$data_raw, "elevation", "dem_china.tif"),
    file.path(DIRS$data_raw, "external", "dem_china.tif")
  )
  for (ep in elev_paths) {
    if (file.exists(ep)) {
      elev_rast <- tryCatch(rast(ep), error = function(e) NULL)
      if (!is.null(elev_rast)) {
        message(sprintf("[03] Loaded DEM from %s", ep))
        break
      }
    }
  }

  # 提取值
  pts <- st_as_sf(centroids, coords = c("centroid_lon", "centroid_lat"), crs = 4326)
  result <- tibble(grid_cell = centroids$grid_cell)

  if (!is.null(bio_stack)) {
    bio_vals <- terra::extract(bio_stack, pts, ID = FALSE)
    for (bv in bio_vars) {
      if (bv %in% names(bio_vals)) {
        result[[bv]] <- bio_vals[[bv]]
      }
    }
    message(sprintf("[03] Extracted %d bio variables", sum(bio_vars %in% names(result))))
  }

  if (!is.null(elev_rast)) {
    elev_vals <- terra::extract(elev_rast, pts, ID = FALSE)
    result$elev_mean <- elev_vals[, 1]
    result$elev_sd <- NA_real_  # 需要区域统计才能得到 sd
    message("[03] Extracted elevation")
  }

  # 添加质心坐标
  result$centroid_lon <- centroids$centroid_lon
  result$centroid_lat <- centroids$centroid_lat

  # 占位：texture_shannon, habitat_diversity_shannon, hfi_mean, landcover_built, landcover_cropland
  placeholder_cols <- c("elev_sd", "texture_shannon", "habitat_diversity_shannon",
                         "hfi_mean", "landcover_built", "landcover_cropland")
  for (pc in placeholder_cols) {
    if (!pc %in% names(result)) result[[pc]] <- NA_real_
  }

  # 如果大量列 NA，从 v2 grid_env 做空间最近邻插值
  v2_env <- safe_read(file.path(DIRS$v2_derived, "grid_environment_dynamic_occupancy.rds"))
  if (!is.null(v2_env)) {
    na_cols <- names(result)[sapply(result, function(x) sum(is.na(x)) > nrow(result) * 0.5)]
    if (length(na_cols) > 0) {
      message(sprintf("[03] Imputing %d NA-heavy columns from v2 grid_env via nearest neighbor",
                      length(na_cols)))
      for (nc in na_cols) {
        if (!nc %in% names(v2_env)) next
        # 对每个 10km 网格，找最近的 v2 100km 网格
        for (i in seq_len(nrow(result))) {
          if (!is.na(result[[nc]][i])) next
          dlat <- abs(v2_env$centroid_lat - result$centroid_lat[i])
          dlon <- abs(v2_env$centroid_lon - result$centroid_lon[i])
          d <- dlat^2 + dlon^2
          nearest <- which.min(d)
          if (length(nearest) > 0 && !is.na(v2_env[[nc]][nearest])) {
            result[[nc]][i] <- v2_env[[nc]][nearest]
          }
        }
      }
    }
  }

  message(sprintf("[03] Raster extraction complete: %d grid cells, %d variables",
                  nrow(result), ncol(result) - 1))
  result
}
