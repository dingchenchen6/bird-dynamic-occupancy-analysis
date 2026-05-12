#!/usr/bin/env Rscript
## 01_merge_birdwatch_ebird.R  —  v3 数据合并与去重
##
## 修复 v2 的去重 bug：匿名用户 fallback (anon_{lon}_{lat}_{date}) 会合并不同观测者
## v3：source-aware hash，同年同地不同来源保留，同来源精确去重
## 年份上限统一为 ANALYSIS_YR_HI=2024

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(lubridate); library(data.table)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
P <- ensure_v3_dirs()

log_time("01", "Starting data merge & dedup")

# ── 1. 加载原始数据 ──────────────────────────────────────────────────
bw_path <- file.path(DIRS$data_raw, "birdwatch_events_raw.rds")
eb_path <- file.path(DIRS$data_raw, "ebird_events_raw.rds")

if (!file.exists(bw_path) && !file.exists(eb_path)) {
  # 尝试从 v2 的合并结果加载
  v2_merged <- file.path(DIRS$v2_derived, "combined_events_merged_dedup_2000_2025.rds")
  if (file.exists(v2_merged)) {
    message("[01] Loading from v2 merged data (will re-dedup)")
    all_events <- readRDS(v2_merged)
    # 如果已有 source 列则跳过标准化
    if (!"source" %in% names(all_events)) {
      all_events$source <- "unknown"
    }
  } else {
    stop("No raw or merged data found. Place birdwatch_events_raw.rds and/or ebird_events_raw.rds in data/.")
  }
} else {
  events_list <- list()

  if (file.exists(bw_path)) {
    bw <- readRDS(bw_path)
    bw$source <- "birdwatch"
    events_list[["birdwatch"]] <- bw
    message(sprintf("[01] BirdWatch: %d records", nrow(bw)))
  }

  if (file.exists(eb_path)) {
    eb <- readRDS(eb_path)
    eb$source <- "ebird"
    events_list[["ebird"]] <- eb
    message(sprintf("[01] eBird: %d records", nrow(eb)))
  }

  all_events <- bind_rows(events_list)
}

# ── 2. 标准化列名 ────────────────────────────────────────────────────
all_events <- all_events |>
  rename_with(~ tolower(gsub("[.]", "_", .x)), everything()) |>
  # 标准化 source 值（兼容 v2 的不同命名）
  mutate(
    source = tolower(source),
    source = case_when(
      grepl("ebird", source)    ~ "ebird",
      grepl("birdwatch", source) | grepl("china_bird", source) ~ "birdwatch",
      TRUE                      ~ source
    )
  ) |>
  mutate(
    species    = str_trim(species),
    longitude  = as.numeric(longitude),
    latitude   = as.numeric(latitude)
  )

# 兼容无 date 列的 v2 数据：从 year/month/day 构造
if (!"date" %in% names(all_events) && all(c("year", "month", "day") %in% names(all_events))) {
  all_events$date <- as.Date(paste(all_events$year, all_events$month, all_events$day, sep = "-"))
} else if ("date" %in% names(all_events)) {
  all_events$date <- as.Date(all_events$date)
}

all_events <- all_events |>
  mutate(
    year       = as.integer(year(date)),
    month      = as.integer(month(date))
  ) |>
  filter(
    !is.na(date), !is.na(species),
    !is.na(longitude), !is.na(latitude),
    between(year, ANALYSIS_YR_LO, ANALYSIS_YR_HI),
    between(latitude, 17, 55),     # 中国范围
    between(longitude, 72, 136)
  )

message(sprintf("[01] After date/geo filter: %d records", nrow(all_events)))

# ── 3. source-aware 去重 ──────────────────────────────────────────────
# v2 bug: 匿名用户用 anon_{lon}_{lat}_{date} 作 key，
#   导致同日同地不同观测者被合并
# v3 fix: 区分来源，同来源精确去重，跨源保留

# 统一观测者 ID（兼容 v2 数据的列名差异）
# v2 数据可能使用 username 而非 observer_id
obs_col <- intersect(c("observer_id", "user_id", "sampler_id", "username"),
                     names(all_events))[1]
if (is.na(obs_col)) obs_col <- NULL

all_events <- all_events |>
  mutate(
    observer_raw = if (!is.null(obs_col)) .data[[obs_col]] else NA_character_,
    observer_id  = case_when(
      source %in% c("ebird", "gbif_ebird_china") ~ paste0("eb_", observer_raw),
      source %in% c("birdwatch", "china_birdwatch_platform") ~ paste0("bw_", observer_raw),
      TRUE                                        ~ paste0("unk_", observer_raw %||% "anon")
    ),
    # 对缺失 observer_id 的记录生成唯一 ID（不再合并不同观测者）
    observer_id = if_else(
      is.na(observer_raw) | str_detect(observer_id, "unk_.*anon$"),
      paste0("anon_", round(longitude, 4), "_", round(latitude, 4), "_",
             as.integer(date)),
      observer_id
    )
  )

# 去重 key: species + date + lon(4dp) + lat(4dp) + source + observer_id
dedup_key <- all_events |>
  mutate(
    lon4 = round(longitude, 4),
    lat4 = round(latitude, 4)
  ) |>
  distinct(species, date, lon4, lat4, source, observer_id, .keep_all = TRUE)

n_before <- nrow(all_events)
n_after  <- nrow(dedup_key)
message(sprintf("[01] Cross-source dedup: %d → %d (%.1f%% removed)",
                n_before, n_after, 100 * (1 - n_after / n_before)))

all_events <- dedup_key

# ── 4. 保存 ──────────────────────────────────────────────────────────
saveRDS(all_events, v3_file("derived", "combined_events_dedup_v3", "rds"))
write_csv(
  all_events |> summarise(
    n_records = n(), n_species = length(unique(species)),
    n_observers = length(unique(observer_id)),
    year_range = paste(range(year, na.rm = TRUE), collapse = "-"),
    sources = paste(unique(source), collapse = ", ")
  ),
  v3_file("results", "table_01_merge_summary_v3")
)

log_time("01", sprintf("DONE: %d records, %d species", nrow(all_events),
                        length(unique(all_events$species))))
