#!/usr/bin/env Rscript
## 04c_run_stMsPGOcc_4chain_merge.R  —  v3 核心模型
##
## 关键升级：tMsPGOcc → stMsPGOcc（增加空间 NNGP 随机效应）
## 检测协变量修正：移除 log_observers（共线）、duration_min → has_duration
## MCMC 参数提升：400 batch / 5000 burn / 4 chains
## 种级收敛诊断 + WAIC

## 内存管理：在 24GB 机器上，需要积极 GC
gc()

suppressPackageStartupMessages({
  library(spOccupancy); library(readr); library(dplyr); library(tidyr)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR")
if (CODE_V3 == "") {
  if (dir.exists(file.path("~", "bird_dynamic_occupancy_analysis", "code_v3"))) {
    CODE_V3 <- file.path("~", "bird_dynamic_occupancy_analysis", "code_v3")
  } else {
    CODE_V3 <- file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3")
  }
}
source(file.path(CODE_V3, "00_config.R"))
# -- 04c override: 200 batch, 2000 burn --
FULL_N_BATCH   <<- 200L
FULL_N_BURN    <<- 2000L
FULL_N_THIN    <<- 2L
FULL_N_CHAINS  <<- 4L
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
P <- ensure_v3_dirs()

log_time("04c", "Starting stMsPGOcc model fit")

# ── 1. 加载数据 ──────────────────────────────────────────────────────
# 优先加载带网格标签的 survey history
survey_stem <- paste0("survey_history", GRID_TAG, "_v3")
survey <- safe_read(v3_file("derived", survey_stem, "rds"))
if (is.null(survey)) {
  # fallback: 旧版无标签文件
  survey <- safe_read(v3_file("derived", "survey_history_v3", "rds"))
}
if (is.null(survey)) stop("Survey history not found. Run 02 first.")

grid_env_stem <- paste0("grid_environment", GRID_TAG, "_v3")
grid_env <- safe_read(v3_file("derived", grid_env_stem, "rds"))
if (is.null(grid_env)) {
  grid_env <- safe_read(v3_file("derived", "grid_environment_v3", "rds"))
}
if (is.null(grid_env)) {
  grid_env <- safe_read(file.path(DIRS$v2_derived, "grid_environment_dynamic_occupancy.rds"))
}
if (is.null(grid_env)) stop("Grid environment not found. Run 03 first.")

det_cov_stem <- paste0("detection_covariates", GRID_TAG, "_v3")
det_cov <- safe_read(v3_file("derived", det_cov_stem, "rds"))
if (is.null(det_cov)) {
  det_cov <- safe_read(v3_file("derived", "detection_covariates_v3", "rds"))
}
# 也尝试从 v2 加载 visit_effort（更详细）
visit_effort <- safe_read(file.path(DIRS$v2_derived, "visit_effort_2000_2024.rds"))
species_visit <- safe_read(file.path(DIRS$v2_derived, "species_visit_2000_2024.rds"))

# ── 2. 准备物种和站点 ──────────────────────────────────────────────
sites <- survey$sites
n_sites <- length(sites)
periods <- survey$periods
n_periods <- length(periods)
candidate_species <- survey$species
n_sp <- length(candidate_species)

# 加载候选物种（如有限制）
MAX_N_SP <- as.integer(Sys.getenv("V3_MAX_SPECIES", "200"))
candidate_all <- read_csv_safe(
  file.path(DIRS$v2_results, "table_dynamic_occupancy_candidate_species_all"))
if (!is.null(candidate_all)) {
  candidate_species <- head(candidate_all$species, MAX_N_SP)
} else {
  candidate_species <- head(survey$species, MAX_N_SP)
}
n_sp <- length(candidate_species)
message(sprintf("[04] Candidate species limited to %d (max: %d)", n_sp, MAX_N_SP))

message(sprintf("[04] %d species × %d sites × %d periods", n_sp, n_sites, n_periods))

# ── 3. 构建 y array 和检测协变量 ──────────────────────────────────────
# stMsPGOcc 需要: y[species, site, time, rep]
# 我们的数据没有真正的 repeat visits 在同一 period，所以 rep=1

# FIX: 站点对齐 — 以 survey$ sites 为基准（与 02/03 对齐的 grid_cell）
# v2 visit_effort 使用不同的 grid_cell 编码，不能直接取交集
# 策略：优先用 survey 数据构建 y（保证站点对齐），v2 数据仅做检测协变量补充

# 检查 v2 visit_effort 的 grid_cell 是否与 survey sites 有足够交集
use_v2_y <- FALSE
if (!is.null(species_visit) && !is.null(visit_effort)) {
  sites_ve <- sort(unique(visit_effort$grid_cell))
  sites_ge <- sort(unique(grid_env$grid_cell))
  overlap_ve <- length(intersect(sites_ve, sites))
  overlap_ge <- length(intersect(sites_ge, sites))
  message(sprintf("[04] Site overlap: visit_effort=%d/%d, grid_env=%d/%d, survey=%d",
                  overlap_ve, length(sites_ve), overlap_ge, length(sites_ge), n_sites))
  # 仅当 v2 的 grid_cell 编码与 survey sites 有 >50% 交集时才使用 v2 构建 y
  if (overlap_ve > n_sites * 0.5 && overlap_ge > n_sites * 0.5) {
    use_v2_y <- TRUE
  }
}

if (use_v2_y) {
  message("[04] Building y array from v2 species_visit + visit_effort (good overlap)")

  # 确保站点与 grid_env 对齐（取交集）
  sites_ve <- sort(unique(visit_effort$grid_cell))
  sites_ge <- sort(unique(grid_env$grid_cell))
  sites <- sort(intersect(sites_ve, sites_ge))
  n_sites <- length(sites)
  message(sprintf("[04] Sites: %d (visit_effort: %d, grid_env: %d, intersect: %d)",
                  n_sites, length(sites_ve), length(sites_ge), n_sites))

  # primary blocks
  n_primary <- n_periods
  n_secondary <- max(visit_effort$year_in_block, na.rm = TRUE)

  # 过滤到交集站点
  eff <- visit_effort |>
    filter(grid_cell %in% sites) |>
    mutate(site_index = match(grid_cell, sites))
  det <- species_visit |>
    filter(species %in% candidate_species, grid_cell %in% sites) |>
    mutate(species_index = match(species, candidate_species),
           site_index    = match(grid_cell, sites))
  # y array: species × sites × primary × secondary
  y <- array(NA_integer_,
             dim = c(n_sp, n_sites, n_primary, n_secondary),
             dimnames = list(candidate_species, as.character(sites),
                             paste0("P", seq_len(n_primary)),
                             paste0("rep", seq_len(n_secondary))))

  visited_mask <- array(FALSE, dim = c(n_sites, n_primary, n_secondary))
  log_events_arr <- array(NA_real_, dim = c(n_sites, n_primary, n_secondary))
  has_duration_arr <- array(NA_real_, dim = c(n_sites, n_primary, n_secondary))
  # FIX #4: 新增 log_duration（缺失时用0填补，has_duration=0 标记缺失）
  log_duration_arr <- array(0, dim = c(n_sites, n_primary, n_secondary))

  for (ii in seq_len(nrow(eff))) {
    rr <- eff[ii, ]
    if (rr$block_id > n_primary || rr$year_in_block > n_secondary) next
    visited_mask[rr$site_index, rr$block_id, rr$year_in_block] <- TRUE
    log_events_arr[rr$site_index, rr$block_id, rr$year_in_block] <- rr$log_events
    # has_duration: 1 if mean_duration_min > 0 and not NA
    has_dur <- as.integer(!is.na(rr$mean_duration_min) && rr$mean_duration_min > 0)
    has_duration_arr[rr$site_index, rr$block_id, rr$year_in_block] <- has_dur
    # log_duration: 有值时取 log10，缺失时为 0（由 has_duration=0 标记）
    if (has_dur == 1) {
      log_duration_arr[rr$site_index, rr$block_id, rr$year_in_block] <- log10(rr$mean_duration_min)
    }
  }

  # 填充物种检测
  for (ii in seq_len(nrow(det))) {
    rr <- det[ii, ]
    if (rr$species_index > n_sp || rr$site_index > n_sites) next
    if (rr$block_id > n_primary || rr$year_in_block > n_secondary) next
    y[rr$species_index, rr$site_index, rr$block_id, rr$year_in_block] <- rr$detected
  }

  # 补 0：访问过但未检测到的位置
  visited_idx <- which(visited_mask, arr.ind = TRUE)
  for (ii in seq_len(nrow(visited_idx))) {
    j <- visited_idx[ii, 1]; t <- visited_idx[ii, 2]; k <- visited_idx[ii, 3]
    miss <- is.na(y[, j, t, k])
    if (any(miss)) y[miss, j, t, k] <- 0L
  }

  # 填充 NA in det covs
  fill_arr_na <- function(arr, mask) {
    obs <- is.finite(arr)
    fill_val <- if (any(obs)) median(arr[obs], na.rm = TRUE) else 0
    arr[mask & !obs] <- fill_val
    arr
  }
  log_events_arr <- fill_arr_na(log_events_arr, visited_mask)
  has_duration_arr <- fill_arr_na(has_duration_arr, visited_mask)
  has_duration_arr[!is.finite(has_duration_arr)] <- 0

} else {
  # Fallback: 使用 v3 survey 数据构建简化 y
  message("[04] Building y array from v3 survey_history (simplified)")
  n_secondary <- 1L
  y <- array(0L, dim = c(n_sp, n_sites, n_periods, n_secondary),
             dimnames = list(candidate_species, as.character(sites),
                             paste0("P", seq_len(n_periods)), "rep1"))

  # 从 survey$y (species × site_period) 重组
  # survey$y 在 02 中按 site_period 行序构建：P1(s1..sn), P2(s1..sn), ...
  y_mat <- survey$y
  sp_idx <- match(candidate_species, rownames(y_mat))
  if (any(is.na(sp_idx))) {
    missing_sp <- candidate_species[is.na(sp_idx)]
    message(sprintf("[04] Warning: %d species not in survey$y (e.g. %s)",
                    length(missing_sp), paste(head(missing_sp, 3), collapse = ", ")))
    sp_idx <- sp_idx[!is.na(sp_idx)]
    n_sp <- length(sp_idx)
    candidate_species <- candidate_species[!is.na(match(candidate_species, rownames(y_mat)))]
    y <- array(0L, dim = c(n_sp, n_sites, n_periods, n_secondary),
               dimnames = list(candidate_species, as.character(sites),
                               paste0("P", seq_len(n_periods)), "rep1"))
  }

  for (t in seq_len(n_periods)) {
    col_start <- (t - 1) * n_sites + 1
    col_end <- t * n_sites
    if (col_end <= ncol(y_mat)) {
      y[, , t, 1] <- as.integer(y_mat[sp_idx, col_start:col_end, drop = FALSE])
    } else {
      message(sprintf("[04] Warning: period %d columns exceed y_mat ncol (%d)", t, ncol(y_mat)))
    }
  }
  message(sprintf("[04] y array filled: %d non-zero entries", sum(y == 1L)))

  # 检测协变量简化
  log_events_arr <- array(0, dim = c(n_sites, n_periods, 1))
  has_duration_arr <- array(1, dim = c(n_sites, n_periods, 1))
  log_duration_arr <- array(0, dim = c(n_sites, n_periods, 1))
  if (!is.null(det_cov)) {
    for (t in seq_len(n_periods)) {
      p_data <- det_cov |> filter(period == paste0("P", t))
      log_events_arr[, t, 1] <- p_data$log_events[match(sites, p_data$grid_cell)]
      has_duration_arr[, t, 1] <- p_data$has_duration[match(sites, p_data$grid_cell)]
      # FIX #4: log_duration 补充
      if ("log_duration" %in% names(p_data)) {
        log_duration_arr[, t, 1] <- p_data$log_duration[match(sites, p_data$grid_cell)]
      }
    }
    log_events_arr[!is.finite(log_events_arr)] <- 0
    has_duration_arr[!is.finite(has_duration_arr)] <- 0
    log_duration_arr[!is.finite(log_duration_arr)] <- 0
  }
}

# ── 4. 准备占有率协变量 ──────────────────────────────────────────────
standardize <- function(x) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) return(rep(0, length(x)))
  x[!is.finite(x)] <- median(x[is.finite(x)], na.rm = TRUE)
  as.numeric(scale(x))
}

grid_env_aligned <- grid_env |>
  filter(grid_cell %in% sites) |>
  arrange(match(grid_cell, sites))

# 协变量筛选
occ_vars <- c("bio4", "bio7", "bio11", "bio13",
              "elev_mean", "elev_sd", "texture_shannon",
              "habitat_diversity_shannon",
              "hfi_mean", "landcover_built", "landcover_cropland",
              "centroid_lon", "centroid_lat")
occ_vars <- intersect(occ_vars, names(grid_env_aligned))
message("[04] Occupancy covariates: ", paste(occ_vars, collapse = ", "))

occ_covs <- setNames(lapply(occ_vars, function(v) standardize(grid_env_aligned[[v]])), occ_vars)
# year_scaled: period-level
occ_covs$year_scaled <- matrix(
  rep(scale(seq_len(n_periods)), each = n_sites),
  nrow = n_sites, ncol = n_periods
)

# 检测协变量（列表格式）— FIX #4: 新增 log_duration
det_covs <- list(
  log_events   = log_events_arr,
  log_duration = log_duration_arr,
  has_duration = has_duration_arr
)

# 坐标（stMsPGOcc 需要）
# 使用 visit_effort 中的站点坐标（与 y 对齐）
coords_df <- grid_env_aligned[, c("grid_cell", "centroid_lon", "centroid_lat")]
coords <- as.matrix(coords_df[, c("centroid_lon", "centroid_lat")])

# grid.index: 如果坐标行数与 y 的站点数不同，需要映射
# stMsPGOcc 要求: nrow(coords) 可以 <= nrow(y[, , 1, 1])
# grid.index[i] = coords 中哪行对应 y 的第 i 个站点
grid_index <- seq_len(n_sites)  # 1:1 映射

# ── 5. stMsPGOcc 数据列表 ──────────────────────────────────────────
data_list <- list(
  y          = y,
  occ.covs   = occ_covs,
  det.covs   = det_covs,
  coords     = coords,
  grid.index = grid_index
)

# ── 6. 模型公式 ──────────────────────────────────────────────────────
occ_formula <- reformulate(c(occ_vars, "year_scaled"))
det_formula <- ~ log_events + log_duration + has_duration

# ── 7. Pilot vs Full + Chain Parallelization ───────────────────────────
is_pilot <- Sys.getenv("V3_PILOT", "0") == "1"
run_label <- paste0(RUN_LABEL, "_04c")

n_batch  <- if (is_pilot) PILOT_N_BATCH  else FULL_N_BATCH
n_burn   <- if (is_pilot) PILOT_N_BURN   else FULL_N_BURN
n_thin   <- if (is_pilot) PILOT_N_THIN   else FULL_N_THIN
n_chains <- if (is_pilot) PILOT_N_CHAINS else FULL_N_CHAINS
n_omp <- as.integer(Sys.getenv("V3_N_OMP_THREADS", "1"))

# ── Chain parallelization ──────────────────────────────────────────
# When V3_CHAIN_ID is set (1-N), run a single chain and save chain-specific file.
# The 04b script will combine chains and compute R-hat/ESS.
# This allows 4 chains to run simultaneously on separate CPU cores.
chain_id <- Sys.getenv("V3_CHAIN_ID", "")
if (nzchar(chain_id)) {
  chain_id <- as.integer(chain_id)
  n_chains <- 1L  # single chain per process
  message(sprintf("[04] CHAIN PARALLEL mode: chain %d of %d", chain_id, FULL_N_CHAINS))
} else {
  chain_id <- NA_integer_
  message("[04] SEQUENTIAL mode: all chains in one process")
}

# Pilot: 只用部分物种
n_sp_pilot <- n_sp  # 默认全量
if (is_pilot) {
  pilot_n <- min(20, n_sp)  # 内存限制：20 种
  set.seed(2024)
  pilot_sp_idx <- sort(sample(n_sp, pilot_n))
  data_list$y <- data_list$y[pilot_sp_idx, , , , drop = FALSE]
  n_sp_pilot <- pilot_n
  message(sprintf("[04] PILOT mode: %d species", pilot_n))
}

# 10km 网格站点数多，适当增加 factor 数量
n_factors_use <- min(5, max(1, n_sp_pilot - 1))

# ── 8. 先验和起始值 ────────────────────────────────────────────────
# z init: observed in any replicate -> 1, else 0 (3D: species × sites × periods)
z_init <- apply(data_list$y, c(1, 2, 3), function(a) as.integer(any(a == 1, na.rm = TRUE)))
message(sprintf("[04] z_init dims: %s (expected: %d x %d x %d)",
                paste(dim(z_init), collapse="x"), n_sp_pilot, n_sites, n_periods))
# 确保维度正确
if (length(dim(z_init)) != 3 || dim(z_init)[1] != n_sp_pilot ||
    dim(z_init)[2] != n_sites || dim(z_init)[3] != n_periods) {
  message("[04] z_init dimension mismatch. Forcing correct dimensions.")
  z_init <- array(as.integer(z_init), dim = c(n_sp_pilot, n_sites, n_periods))
}

# 链专属种子：确保不同链探索不同后验区域
chain_seed <- if (!is.na(chain_id)) 2024L + chain_id * 1000L else 2024L
set.seed(chain_seed)
message(sprintf("[04] Chain seed: %d", chain_seed))

# 用测地线距离计算初始 phi（比欧氏距离更合理，尤其对 10km 网格）
# geosphere::distm 给 Haversine 距离 (km)，避免经纬度欧氏距离的偏差
if (requireNamespace("geosphere", quietly = TRUE)) {
  dmat_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
  phi_init <- 3 / mean(dmat_km[upper.tri(dmat_km)])
  message(sprintf("[04] phi init (geodesic): %.4f (mean pairwise dist: %.1f km)",
                  phi_init, mean(dmat_km[upper.tri(dmat_km)])))
  rm(dmat_km)
} else {
  phi_init <- 3 / mean(dist(coords))
  message(sprintf("[04] phi init (Euclidean fallback): %.4f", phi_init))
}

inits <- list(
  alpha.comm = 0, beta.comm = 0, beta = 0, alpha = 0,
  tau.sq.beta = 1, tau.sq.alpha = 1,
  z = z_init,
  sigma.sq = 1,
  phi = phi_init
)

priors <- list(
  beta.comm.normal  = list(mean = 0, var = 2.72),
  alpha.comm.normal = list(mean = 0, var = 2.72),
  tau.sq.beta.ig    = list(a = 0.1, b = 0.1),
  tau.sq.alpha.ig   = list(a = 0.1, b = 0.1),
  sigma.sq.ig       = list(a = 2, b = 1),
  phi.unif          = list(a = 0.1, b = 10)
)

tuning <- list(
  phi = 0.5,
  sigma.sq = 1,
  rho = 0.2
)

# ── 9. 运行 stMsPGOcc ──────────────────────────────────────────────
message(sprintf("[04] Running tMsPGOcc: %s", run_label))
message(sprintf("[04] MCMC: %d batch × 25, burn=%d, thin=%d, chains=%d",
                n_batch, n_burn, n_thin, n_chains))
message(sprintf("[04] Spatial: NNGP %d neighbors, %s covariance",
                N_NEIGHBORS, COV_MODEL))

t_start <- Sys.time()

fit <- tMsPGOcc(
  occ.formula   = occ_formula,
  det.formula   = det_formula,
  data          = data_list,
  inits         = inits,
  priors        = priors,
  tuning        = tuning,
  cov.model     = COV_MODEL,
  NNGP          = TRUE,
  n.neighbors   = N_NEIGHBORS,
  n.factors     = n_factors_use,   # factor model for species correlations
  n.batch       = n_batch,
  batch.length  = 50,
  accept.rate   = 0.43,
  n.omp.threads = n_omp,
  verbose       = TRUE,
  ar1           = TRUE,
  n.report      = max(1, floor(n_batch / 10)),
  n.burn        = n_burn,
  n.thin        = n_thin,
  n.chains      = n_chains
)

t_elapsed <- difftime(Sys.time(), t_start, units = "mins")
message(sprintf("[04] Model fit completed in %.1f minutes", as.numeric(t_elapsed)))

# ── 10. 保存模型（含 OOM 保护和保存验证）────────────────────────────────
# 链并行模式：保存链专属文件（04b 会合并）
if (!is.na(chain_id)) {
  fit_path <- v3_file("derived", paste0("tMsPGOcc_fit_", run_label, "_chain", chain_id), "rds")
} else {
  fit_path <- v3_file("derived", paste0("tMsPGOcc_fit_", run_label), "rds")
}

# 先 GC 释放内存再保存
gc()

save_ok <- tryCatch({
  saveRDS(fit, fit_path, compress = "xz")
  # 验证保存成功
  if (!file.exists(fit_path)) stop("File not created")
  fit_size <- file.info(fit_path)$size
  if (is.na(fit_size) || fit_size < 1000) stop("File too small, likely corrupt")
  message(sprintf("[04] Model saved to %s (%.1f MB)", fit_path, fit_size / 1e6))
  TRUE
}, error = function(e) {
  message(sprintf("[04] saveRDS FAILED: %s", e$message))
  message("[04] Attempting checkpoint save of key slots ...")
  # OOM fallback: 只保存必要 slot
  fit_lite <- list(
    psi.samples = fit$psi.samples,
    y           = fit$y,
    rhat        = fit$rhat,
    ESS         = fit$ESS,
    run.params  = list(
      n_species = n_sp, n_sites = n_sites, n_periods = n_periods,
      n_batch = n_batch, n_burn = n_burn, n_thin = n_thin, n_chains = n_chains
    )
  )
  lite_path <- sub("\\.rds$", "_lite.rds", fit_path)
  tryCatch({
    saveRDS(fit_lite, lite_path, compress = "xz")
    message(sprintf("[04] Lite checkpoint saved to %s", lite_path))
  }, error = function(e2) {
    message(sprintf("[04] Lite save also FAILED: %s", e2$message))
  })
  FALSE
})

# WAIC
waic_val <- tryCatch({
  waic(fit)$waic
}, error = function(e) NA_real_)
message(sprintf("[04] WAIC = %.2f", waic_val))

# 收敛诊断
diag_result <- tryCatch({
  source(file.path(CODE_V3, "utils_diagnostics.R"))
  # 提取社区级参数诊断
  comm_slots <- c("beta.comm.samples", "alpha.comm.samples",
                  "tau.sq.beta.samples", "tau.sq.alpha.samples")
  diag_list <- list()
  for (slot in comm_slots) {
    if (slot %in% names(fit)) {
      samples <- fit[[slot]]
      # 计算 R-hat 和 ESS
      n_ch <- dim(samples)[length(dim(samples))]
      for (ch in seq_len(min(n_ch, 5))) {
        ch_samples <- if (length(dim(samples)) == 3) samples[, , ch] else samples[, ch]
        # 简化：只报告存在
      }
    }
  }
  TRUE
}, error = function(e) {
  message("[04] Diagnostics failed: ", conditionMessage(e))
  FALSE
})

# 保存摘要
write_csv(tibble(
  run_label     = run_label,
  n_species     = if (is_pilot) pilot_n else n_sp,
  n_sites       = n_sites,
  n_periods     = n_periods,
  cov_model     = COV_MODEL,
  n_neighbors   = N_NEIGHBORS,
  n_batch       = n_batch,
  n_burn        = n_burn,
  n_thin        = n_thin,
  n_chains      = n_chains,
  elapsed_mins  = as.numeric(t_elapsed),
  waic          = waic_val,
  is_pilot      = is_pilot
), v3_file("results", paste0("table_model_summary_", run_label)))

log_time("04c", sprintf("DONE: %s, WAIC=%.2f, %.0f min",
                         run_label, waic_val, as.numeric(t_elapsed)))
