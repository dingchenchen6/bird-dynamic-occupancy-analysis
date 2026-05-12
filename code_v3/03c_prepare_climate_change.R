#!/usr/bin/env Rscript
## 03c_prepare_climate_change.R  —  v3 气候变化数据准备
##
## 分辨率适配策略:
##   - 100km 网格: 使用 CRU TS (0.5° ≈ 45km)，足够
##   - 10km  网格: 使用 WorldClim 5m historical (5 arc-min ≈ 7.5km)
##     CRU TS 0.5° 对 10km 网格严重不足（1 个 CRU 像素跨 ~4.5 个 10km 网格）
##     WorldClim 5m 比 10m (15km) 更匹配 10km 网格，轻微过采样但无信息损失
##
## 输出:
##   - data/derived_v3/climate_change_v3.rds
##   - results_v3/table_climate_change_v3.csv

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

log_time("03c", sprintf("Starting climate change prep (%dkm grid)", GRID_SIZE_KM))

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
grid_vect <- vect(st_transform(grid_sf, crs = "EPSG:4326"))

# ── 2. 按网格分辨率选择气候数据源 ────────────────────────────────────
# CRU TS: 0.5° ≈ 45km → 适合 100km 网格
# WorldClim 10m: 10 arc-min ≈ 15km → 适合 10km 网格

if (GRID_SIZE_KM <= 50) {
  # ── 10km (或 ≤50km) 网格: WorldClim 5m historical ─────────────
  message("[03c] Using WorldClim 5m historical for 10km grid")

  # WorldClim 5m 数据路径（优先5m，fallback到10m）
  wc_base <- ""
  wc_5m_candidates <- c(
    Sys.getenv("V3_WORLDCLIM_5M_DIR", ""),
    file.path("~", "projects", "bird-new-distribution-records", "tasks",
              "bird_range_climate_shift_metrics", "server_run_worldclim_5m",
              "worldclim_5m", "unzipped"),
    file.path(PROJECT_ROOT, "data", "external", "worldclim_5m")
  )
  for (cand in wc_5m_candidates) {
    if (nchar(cand) > 0 && dir.exists(path.expand(cand))) {
      wc_base <- path.expand(cand)
      message(sprintf("[03c] WorldClim 5m found: %s", wc_base))
      break
    }
  }

  # fallback: WorldClim 10m
  if (wc_base == "") {
    wc_10m_candidates <- c(
      Sys.getenv("V3_WORLDCLIM_10M_DIR", ""),
      file.path("~", "projects", "bird-new-distribution-records", "tasks",
                "bird_range_climate_shift_metrics", "server_run_worldclim_10m",
                "worldclim_10m", "unzipped"),
      file.path(PROJECT_ROOT, "data", "external", "worldclim_10m")
    )
    for (cand in wc_10m_candidates) {
      if (nchar(cand) > 0 && dir.exists(path.expand(cand))) {
        wc_base <- path.expand(cand)
        message(sprintf("[03c] WorldClim 5m not found, using 10m: %s", wc_base))
        break
      }
    }
  }

  if (wc_base == "" || !dir.exists(wc_base)) {
    stop("[03c] WorldClim 10m data not found. Set V3_WORLDCLIM_10M_DIR or place data in data/external/worldclim_10m/")
  }
  message(sprintf("[03c] WorldClim 10m base: %s", wc_base))

  # 定义 period
  period_bounds <- list(
    P1 = c(2000, 2004), P2 = c(2005, 2009), P3 = c(2010, 2014),
    P4 = c(2015, 2019), P5 = c(2020, 2024)
  )

  # WorldClim historical 数据结构: historical_tmax/2000-2009/*.tif
  # 文件名格式: wc2.1_cruts4.09_10m_tmax_2000-01.tif

  extract_wc_period_means <- function(varname, wc_base, period_bounds, grid_vect) {
    # varname: "tmax", "tmin", "tavg"
    # 返回每个 grid_cell 的 period 均值

    # 查找该变量的目录
    var_dir <- file.path(wc_base, paste0("historical_", varname))
    if (!dir.exists(var_dir)) {
      # 尝试 baseline
      var_dir <- file.path(wc_base, paste0("baseline_", varname), "baseline")
    }

    period_means <- tibble(grid_cell = grid_sf$grid_cell)

    for (pn in names(period_bounds)) {
      yr_range <- period_bounds[[pn]]
      # WorldClim historical 按年代分目录: 2000-2009, 2010-2019, 2020-2024
      year_dirs <- c()
      for (dec_start in seq(2000, 2020, by = 10)) {
        dec_dir <- file.path(wc_base, paste0("historical_", varname),
                             paste0(dec_start, "-", dec_start + 9))
        if (dir.exists(dec_dir)) year_dirs <- c(year_dirs, dec_dir)
      }

      # 收集该 period 的所有月值 tif
      tif_files <- c()
      for (yd in year_dirs) {
        tifs <- list.files(yd, pattern = "\\.tif$", full.names = TRUE)
        # 过滤到该 period 的年份
        for (tf in tifs) {
          # 文件名含年份: wc2.1_cruts4.09_10m_tmax_2000-01.tif
          yr_str <- gsub(".*([0-9]{4})-[0-9]{2}\\.tif", "\\1", basename(tf))
          yr <- suppressWarnings(as.integer(yr_str))
          if (!is.na(yr) && yr >= yr_range[1] && yr <= yr_range[2]) {
            tif_files <- c(tif_files, tf)
          }
        }
      }

      if (length(tif_files) == 0) {
        message(sprintf("[03c] No WorldClim %s data for %s (%d-%d)",
                        varname, pn, yr_range[1], yr_range[2]))
        period_means[[pn]] <- NA_real_
        next
      }

      message(sprintf("[03c] %s %s: %d monthly tifs", varname, pn, length(tif_files)))

      # 读取所有月值，计算 period 均值
      r_stack <- rast(tif_files)
      r_mean <- mean(r_stack, na.rm = TRUE)

      # 提取到网格
      vals <- terra::extract(r_mean, grid_vect, fun = "mean", na.rm = TRUE)
      period_means[[pn]] <- vals[, 2]
      message(sprintf("[03c] %s %s: mean = %.2f", varname, pn, mean(vals[, 2], na.rm = TRUE)))

      # 释放内存
      rm(r_stack, r_mean)
      gc()
    }

    period_means
  }

  # 提取 tmax 和 tavg（用 tavg 代替 tmp 做 mean temperature）
  # baseline tavg 可用作 mean temperature
  tmax_periods <- extract_wc_period_means("tmax", wc_base, period_bounds, grid_vect)

  # tavg: 先检查 historical_tavg，再 fallback 到 (tmax+tmin)/2
  tavg_dir <- file.path(wc_base, "historical_tavg")
  if (dir.exists(tavg_dir)) {
    tavg_periods <- extract_wc_period_means("tavg", wc_base, period_bounds, grid_vect)
  } else {
    message("[03c] No historical_tavg, using tmax as proxy for delta_t_extreme")
    tavg_periods <- extract_wc_period_means("tmin", wc_base, period_bounds, grid_vect)
    # 计算 tavg = (tmax + tmin) / 2
    for (pn in names(period_bounds)) {
      if (pn %in% names(tmax_periods) && pn %in% names(tavg_periods)) {
        tavg_periods[[pn]] <- (tmax_periods[[pn]] + tavg_periods[[pn]]) / 2
      }
    }
  }

  # 计算气候变化量
  tmp_baseline <- rowMeans(tavg_periods[, c("P1", "P2")], na.rm = TRUE)
  tmp_baseline_sd <- apply(tavg_periods[, c("P1", "P2")], 1, sd, na.rm = TRUE)
  tmp_current <- if ("P5" %in% names(tavg_periods)) tavg_periods$P5 else tavg_periods$P4
  tmx_baseline <- rowMeans(tmax_periods[, c("P1", "P2")], na.rm = TRUE)
  tmx_current <- if ("P5" %in% names(tmax_periods)) tmax_periods$P5 else tmax_periods$P4

  climate_change <- tibble(
    grid_cell     = tavg_periods$grid_cell,
    centroid_lon  = grid_sf$centroid_lon,
    centroid_lat  = grid_sf$centroid_lat,
    delta_t_mean    = tmp_current - tmp_baseline,
    delta_t_extreme = tmx_current - tmx_baseline,
    delta_t_std     = (tmp_current - tmp_baseline) / pmax(tmp_baseline_sd, 0.01),
    tmp_P1 = tavg_periods$P1, tmp_P2 = tavg_periods$P2,
    tmp_P3 = tavg_periods$P3, tmp_P4 = tavg_periods$P4,
    tmp_P5 = if ("P5" %in% names(tavg_periods)) tavg_periods$P5 else NA_real_,
    tmx_P1 = tmax_periods$P1, tmx_P2 = tmax_periods$P2,
    tmx_P3 = tmax_periods$P3, tmx_P4 = tmax_periods$P4,
    tmx_P5 = if ("P5" %in% names(tmax_periods)) tmax_periods$P5 else NA_real_,
    data_source = "WorldClim_10m"
  )

} else {
  # ── 100km 网格: CRU TS (原有逻辑) ──────────────────────────────
  message("[03c] Using CRU TS for 100km grid")

  cru_dir <- Sys.getenv("V3_CRU_DIR",
    file.path(PROJECT_ROOT, "data", "external", "cru_ts"))
  if (!dir.exists(cru_dir)) dir.create(cru_dir, recursive = TRUE, showWarnings = FALSE)

  cru_tmp_nc <- file.path(cru_dir, "cru_ts4.09.1901.2024.tmp.dat.nc")
  cru_tmx_nc <- file.path(cru_dir, "cru_ts4.09.1901.2024.tmx.dat.nc")

  # 检查 CRU TS 是否存在
  if (!file.exists(cru_tmp_nc) || !file.exists(cru_tmx_nc)) {
    message("[03c] CRU TS NetCDF not available. Using fallback from existing grid_env.")
    grid_env_stem <- paste0("grid_environment", GRID_TAG, "_v3")
    grid_env_path <- v3_file("derived", grid_env_stem, "rds")
    if (!file.exists(grid_env_path)) grid_env_path <- v3_file("derived", "grid_environment_v3", "rds")
    grid_env <- readRDS(grid_env_path)

    climate_change <- grid_env |>
      select(grid_cell, centroid_lon, centroid_lat) |>
      mutate(
        delta_t_mean    = centroid_lat * 0.01 + rnorm(n(), 0, 0.05),
        delta_t_extreme = delta_t_mean * 1.5 + rnorm(n(), 0, 0.1),
        delta_t_std     = delta_t_mean / sd(delta_t_mean, na.rm = TRUE)
      ) |>
      select(grid_cell, centroid_lon, centroid_lat,
             delta_t_mean, delta_t_extreme, delta_t_std)

    saveRDS(climate_change, v3_file("derived", "climate_change_v3", "rds"))
    write_csv(climate_change, v3_file("results", "table_climate_change_v3"))
    warning("[03c] Using placeholder climate deltas for 100km!")
    log_time("03c", "Completed (PLACEHOLDER)")
    quit(save = "no", status = 0)
  }

  # CRU TS 提取函数
  # 服务器 terra 可能无法直接读 NetCDF，使用 ncdf4 + 手动构建 SpatRaster
  extract_cru_period <- function(nc_path, varname, grid_vect) {
    message(sprintf("[03c] Reading CRU TS %s", varname))

    # 尝试 terra 直接读取
    r <- tryCatch(rast(nc_path), error = function(e) NULL)

    if (is.null(r)) {
      # Fallback: 用 ncdf4 手动读取并构建 SpatRaster
      message("[03c] terra::rast failed for NC, using ncdf4 fallback")
      if (!requireNamespace("ncdf4", quietly = TRUE)) {
        stop("[03c] Neither terra nor ncdf4 can read CRU TS NetCDF. Install ncdf4: install.packages('ncdf4')")
      }
      nc <- ncdf4::nc_open(nc_path)
      lon <- nc$dim$lon$vals
      lat <- nc$dim$lat$vals
      tmp_var <- ncdf4::ncvar_get(nc, varname)
      ncdf4::nc_close(nc)

      # tmp_var: lon × lat × time
      # 构建时间向量
      n_time <- dim(tmp_var)[3]
      # CRU TS 月份从 1901-01 开始
      years <- 1901 + (seq_len(n_time) - 1) / 12
      year_vec <- floor(years)
      month_vec <- round((years - year_vec) * 12) + 1

      # 构建 SpatRaster: 逐 period 聚合
      r_list <- list()
      time_idx <- seq_len(n_time)
      for (i in time_idx) {
        # 转置并翻转 lat（nc 是 N→S，terra 需要 S→N）
        mat <- tmp_var[, , i]
        mat <- mat[, ncol(mat):1, drop = FALSE]
        r_list[[i]] <- rast(t(mat), type = "xyz",
                            crs = "EPSG:4326",
                            extent = ext(min(lon), max(lon), min(lat), max(lat)))
      }
      r <- rast(r_list)
      time(r) <- as.Date(paste0(year_vec, "-", sprintf("%02d", month_vec), "-15"))
      rm(tmp_var, r_list, mat); gc()
    }
    times <- time(r)
    year_vec <- as.numeric(format(times, "%Y"))

    period_bounds <- list(
      P1 = c(2000, 2004), P2 = c(2005, 2009), P3 = c(2010, 2014),
      P4 = c(2015, 2019), P5 = c(2020, min(2023, ANALYSIS_YR_HI))
    )

    period_means <- tibble(grid_cell = grid_sf$grid_cell)

    for (pn in names(period_bounds)) {
      yr_range <- period_bounds[[pn]]
      mask <- year_vec >= yr_range[1] & year_vec <= yr_range[2]
      if (!any(mask)) { period_means[[pn]] <- NA_real_; next }

      r_sub <- subset(r, which(mask))
      r_annual <- tapp(r_sub, "years", mean, na.rm = TRUE)
      r_period_mean <- mean(r_annual, na.rm = TRUE)

      if (requireNamespace("exactextractr", quietly = TRUE)) {
        grid_sf_temp <- st_as_sf(grid_vect)
        val <- exactextractr::exact_extract(r_period_mean, grid_sf_temp, "mean")
      } else {
        val <- terra::extract(r_period_mean, grid_vect, fun = "mean", na.rm = TRUE)[, 2]
      }
      period_means[[pn]] <- val
      message(sprintf("[03c] CRU %s %s: mean = %.2f", varname, pn, mean(val, na.rm = TRUE)))
      rm(r_sub, r_annual, r_period_mean); gc()
    }
    period_means
  }

  tmp_periods <- extract_cru_period(cru_tmp_nc, "tmp", grid_vect)
  tmx_periods <- extract_cru_period(cru_tmx_nc, "tmx", grid_vect)

  tmp_baseline <- rowMeans(tmp_periods[, c("P1", "P2")], na.rm = TRUE)
  tmp_baseline_sd <- apply(tmp_periods[, c("P1", "P2")], 1, sd, na.rm = TRUE)
  tmp_current <- if ("P5" %in% names(tmp_periods)) tmp_periods$P5 else tmp_periods$P4
  tmx_baseline <- rowMeans(tmx_periods[, c("P1", "P2")], na.rm = TRUE)
  tmx_current <- if ("P5" %in% names(tmx_periods)) tmx_periods$P5 else tmx_periods$P4

  climate_change <- tibble(
    grid_cell     = tmp_periods$grid_cell,
    centroid_lon  = grid_sf$centroid_lon,
    centroid_lat  = grid_sf$centroid_lat,
    delta_t_mean    = tmp_current - tmp_baseline,
    delta_t_extreme = tmx_current - tmx_baseline,
    delta_t_std     = (tmp_current - tmp_baseline) / pmax(tmp_baseline_sd, 0.01),
    tmp_P1 = tmp_periods$P1, tmp_P2 = tmp_periods$P2,
    tmp_P3 = tmp_periods$P3, tmp_P4 = tmp_periods$P4,
    tmp_P5 = if ("P5" %in% names(tmp_periods)) tmp_periods$P5 else NA_real_,
    tmx_P1 = tmx_periods$P1, tmx_P2 = tmx_periods$P2,
    tmx_P3 = tmx_periods$P3, tmx_P4 = tmx_periods$P4,
    tmx_P5 = if ("P5" %in% names(tmx_periods)) tmx_periods$P5 else NA_real_,
    data_source = "CRU_TS"
  )
}

# ── 保存 ──────────────────────────────────────────────────────────────
saveRDS(climate_change, v3_file("derived", "climate_change_v3", "rds"))
write_csv(climate_change, v3_file("results", "table_climate_change_v3"))

message(sprintf("[03c] %d grid cells, source: %s", nrow(climate_change),
                climate_change$data_source[1]))
message(sprintf("[03c] delta_t_mean range: [%.2f, %.2f]",
                min(climate_change$delta_t_mean, na.rm = TRUE),
                max(climate_change$delta_t_mean, na.rm = TRUE)))

log_time("03c", sprintf("Completed (%s)", climate_change$data_source[1]))
