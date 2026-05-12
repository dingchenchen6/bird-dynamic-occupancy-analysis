#!/usr/bin/env Rscript
## 02_build_survey_history.R  —  v3 构建调查历史
##
## 关键修复：繁殖季过滤（4-8月，可配置）
## v2 没有月份过滤，包含越冬/迁徙鸟混淆占有信号

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(sf)
  library(data.table)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_spatial.R"))
P <- ensure_v3_dirs()

log_time("02", "Starting survey history build")

# ── 1. 加载去重后的事件数据 ───────────────────────────────────────────
events_path <- v3_file("derived", "combined_events_dedup_v3", "rds")
if (!file.exists(events_path)) {
  # fallback: 尝试 v2 数据
  events_path <- file.path(DIRS$v2_derived, "combined_events_merged_dedup_2000_2025.rds")
}
events <- safe_read(events_path)
if (is.null(events)) stop("No dedup events found. Run 01 first.")

# ── 2. 繁殖季过滤 ────────────────────────────────────────────────────
if (!"month" %in% names(events)) {
  events$month <- as.integer(format(as.Date(events$date), "%m"))
}

events_filtered <- events |>
  filter(month %in% BREEDING_MONTHS)

message(sprintf("[02] Breeding season filter (months %s): %d → %d (%.1f%% retained)",
                paste(BREEDING_MONTHS, collapse = "-"),
                nrow(events), nrow(events_filtered),
                100 * nrow(events_filtered) / nrow(events)))

# FIX #2: 繁殖季过滤审计 — 按 period 统计数据量下降比例
# 如果早期主期损失 >50% 的数据，需要考虑全年建模
# periods 和 n_periods 需在此处先定义（后面也用）
periods <- seq(ANALYSIS_YR_LO, ANALYSIS_YR_HI, by = PERIOD_LENGTH)
n_periods <- length(periods)

events$period_raw <- cut(events$year,
  breaks = c(periods - 0.5, ANALYSIS_YR_HI + 0.5),
  labels = paste0("P", seq_len(n_periods)),
  right = TRUE)
breeding_audit <- events |>
  filter(!is.na(period_raw)) |>
  group_by(period_raw) |>
  summarise(n_total = n(), .groups = "drop")

events_filtered_pre <- events |>
  filter(!is.na(period_raw), month %in% BREEDING_MONTHS) |>
  group_by(period_raw) |>
  summarise(n_breeding = n(), .groups = "drop")

breeding_audit <- breeding_audit |>
  left_join(events_filtered_pre, by = "period_raw") |>
  mutate(pct_retained = 100 * n_breeding / n_total)

message("[02] Breeding filter audit by period:")
for (i in seq_len(nrow(breeding_audit))) {
  r <- breeding_audit[i, ]
  flag <- if (r$pct_retained < 50) " ⚠️ <50%!" else ""
  message(sprintf("  %s: %d → %d (%.1f%% retained)%s",
                  r$period_raw, r$n_total, r$n_breeding, r$pct_retained, flag))
}
write_csv(breeding_audit, v3_file("results", "table_breeding_filter_audit"))
rm(breeding_audit, events_filtered_pre)

# ── 3. 构建网格 × 物种 × period 的调查历史 ────────────────────────────
# 5 年一个 primary period
periods <- seq(ANALYSIS_YR_LO, ANALYSIS_YR_HI, by = PERIOD_LENGTH)
n_periods <- length(periods)

events_filtered <- events_filtered |>
  mutate(
    period = cut(year, breaks = c(periods - 0.5, ANALYSIS_YR_HI + 0.5),
                  labels = paste0("P", seq_len(n_periods)),
                  right = TRUE)
  ) |>
  filter(!is.na(period))

# 加载/构建网格
# 仅在 100km 且有 v2 数据时加载旧网格；否则从头构建
grid_rds_path <- v3_file("derived", paste0("china_grid_", GRID_SIZE_KM, "km_v3"), "rds")
grid_sf <- NULL

if (GRID_SIZE_KM == 100L) {
  # 尝试加载 v2 100km 网格
  grid_sf <- safe_read(file.path(DIRS$v2_derived, "china_grid_100km_v2.rds"))
}

if (is.null(grid_sf) && file.exists(grid_rds_path)) {
  # 尝试加载之前构建的同分辨率网格
  grid_sf <- safe_read(grid_rds_path)
  if (!is.null(grid_sf)) message(sprintf("[02] Loaded existing %dkm grid", GRID_SIZE_KM))
}

if (is.null(grid_sf)) {
  message(sprintf("[02] Building %dkm grid from scratch", GRID_SIZE_KM))
  # 从数据点构建网格
  events_sf <- st_as_sf(events_filtered,
                         coords = c("longitude", "latitude"),
                         crs = st_crs(4326))
  # Albers 投影
  events_albers <- project_china_albers(events_sf)

  # 按配置的网格大小构建
  bbox <- st_bbox(events_albers)
  grid_sf <- st_make_grid(events_albers, cellsize = GRID_SIZE_KM * 1000,
                            square = TRUE) |>
    st_sf() |>
    mutate(grid_cell = row_number())
  grid_sf <- st_transform(grid_sf, st_crs(4326))

  # 保存网格供后续使用
  saveRDS(grid_sf, grid_rds_path)
  message(sprintf("[02] Grid saved to %s (%d cells)", grid_rds_path, nrow(grid_sf)))
}

# 空间匹配
events_sf <- st_as_sf(events_filtered,
                       coords = c("longitude", "latitude"),
                       crs = st_crs(4326))
sf_use_s2(FALSE)
grid_match <- st_join(events_sf, grid_sf, join = st_intersects)

# 移除未匹配的
grid_match <- grid_match |>
  filter(!is.na(grid_cell)) |>
  st_drop_geometry()

# ── 4. 构建物种列表（按频率筛选）───────────────────────────────────────
species_counts <- grid_match |>
  count(species, sort = TRUE) |>
  filter(n >= MIN_VISITS)

candidate_species <- species_counts$species
message(sprintf("[02] %d candidate species (≥ %d records)", length(candidate_species), MIN_VISITS))

# ── 5. 构建 y (detection/nondetection) 矩阵 ────────────────────────────
# 行 = 物种，列 = 网格 × period 组合
sites <- sort(unique(grid_match$grid_cell))
site_period <- expand.grid(grid_cell = sites, period = paste0("P", seq_len(n_periods)),
                            stringsAsFactors = FALSE)
site_period$site_period_id <- seq_len(nrow(site_period))

# 访问矩阵
visit_matrix <- grid_match |>
  count(grid_cell, period, name = "visits") |>
  right_join(site_period, by = c("grid_cell", "period")) |>
  mutate(visits = replace_na(visits, 0L))

# 检测矩阵
y_list <- list()
for (sp in candidate_species) {
  sp_dat <- grid_match |>
    filter(species == sp) |>
    distinct(grid_cell, period, .keep_all = FALSE) |>
    mutate(detected = 1L)
  sp_full <- site_period |>
    left_join(sp_dat, by = c("grid_cell", "period")) |>
    mutate(detected = replace_na(detected, 0L))
  y_list[[sp]] <- sp_full$detected
}
y_mat <- do.call(rbind, y_list)
rownames(y_mat) <- candidate_species

# 访问次数矩阵（spOccupancy 格式）
occ_visits <- visit_matrix |>
  pivot_wider(names_from = period, values_from = visits,
               values_fill = 0) |>
  arrange(match(grid_cell, sites))

# ── 6. 保存 ──────────────────────────────────────────────────────────
survey_stem <- paste0("survey_history", GRID_TAG, "_v3")
saveRDS(list(
  y              = y_mat,
  visits         = occ_visits,
  sites          = sites,
  species        = candidate_species,
  periods        = periods,
  site_period    = site_period,
  breeding_months = BREEDING_MONTHS,
  grid_size_km   = GRID_SIZE_KM
), v3_file("derived", survey_stem, "rds"))

write_csv(tibble(
  n_species = length(candidate_species),
  n_sites   = length(sites),
  n_periods = n_periods,
  breeding_months = paste(BREEDING_MONTHS, collapse = ","),
  n_events_total  = nrow(events),
  n_events_breeding = nrow(events_filtered),
  grid_size_km = GRID_SIZE_KM
), v3_file("results", paste0("table_02_survey_summary", GRID_TAG, "_v3")))

log_time("02", sprintf("DONE: %d sp × %d sites × %d periods",
                        length(candidate_species), length(sites), n_periods))
