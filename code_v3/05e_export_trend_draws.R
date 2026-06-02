#!/usr/bin/env Rscript
## 05e_export_trend_draws.R  —  从 metric_arrays_checkpoint 直接生成 trend_draws_<metric>_*.rds
##
## 背景:05_postprocess_diversity_extended.R 在内存里算了 trend_draws 但没 saveRDS,
## 导致 09 [C2] RF 重要性必须降级为 single-draw(无后验传播)。
## 本脚本读取已有 metric_arrays_checkpoint,直接重算并落盘 8 个 metric × Theil-Sen 趋势抽样。
## 输出格式与 09/utils_importance::rf_importance_draws 期望一致:grid × draws 矩阵。
##
## 输入:data/derived_v3/metric_arrays_checkpoint_<run_label>_extended.rds
## 输出:data/derived_v3/trend_draws_<metric>_<run_label>_extended.rds × 8

suppressPackageStartupMessages({
  library(parallel)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_diversity.R"))  # theil_sen_slope

P <- ensure_v3_dirs()

is_pilot <- Sys.getenv("V3_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL
run_label_ext <- paste0(run_label, "_extended")

log_time("05e", "Export per-metric Theil-Sen trend_draws from checkpoint")

# ── 1. 加载 metric checkpoint ───────────────────────────────────────
ck_path <- v3_file("derived", paste0("metric_arrays_checkpoint_", run_label_ext), "rds")
if (!file.exists(ck_path)) {
  ck_path_alt <- v3_file("derived", paste0("metric_arrays_", run_label_ext), "rds")
  if (file.exists(ck_path_alt)) {
    ck_path <- ck_path_alt
  } else {
    stop(sprintf("metric_arrays_checkpoint not found: %s", ck_path))
  }
}
message(sprintf("[05e] Loading: %s", ck_path))
ck <- readRDS(ck_path)

metric_arr <- ck$metric_arr
sites <- ck$sites
metrics <- ck$metrics
n_periods <- ck$n_periods
n_draws <- dim(metric_arr)[1]
n_sites <- length(sites)

message(sprintf("[05e] metric_arr dim: %d draws x %d sites x %d periods x %d metrics",
  n_draws, n_sites, n_periods, length(metrics)))
message(sprintf("[05e] available metrics: %s", paste(metrics, collapse = ", ")))

# ── 2. 待导出的 metric 集合 ─────────────────────────────────────────
metrics_export <- c("corrected_richness", "shannon", "inv_simpson",
                    "fric_prob", "fdis_prob", "feve_fund",
                    "pd_prob_mctavish", "mpd_prob_mctavish")
metrics_export <- intersect(metrics_export, metrics)
message(sprintf("[05e] exporting %d metrics: %s",
  length(metrics_export), paste(metrics_export, collapse = ", ")))

# ── 3. 计算 Theil-Sen 趋势抽样(draw × site)然后转置 ────────────────
mc_cores <- if (exists("MC_CORES")) MC_CORES else max(1L, parallel::detectCores() - 2L)
mc_cores <- min(mc_cores, 16L)  # 防止超量 fork

for (m in metrics_export) {
  t_start <- Sys.time()
  message(sprintf("[05e] -> metric: %s (cores=%d)", m, mc_cores))

  # 取该 metric 的 [draw, site, period] 子立方
  arr_m <- metric_arr[, , , m, drop = FALSE]
  # arr_m dim:[n_draws, n_sites, n_periods, 1]

  # 并行:每个 draw 一次,内部对全部 site 做 Theil-Sen
  td_list <- mclapply(seq_len(n_draws), function(d) {
    v <- numeric(n_sites)
    for (s in seq_len(n_sites)) {
      vals <- arr_m[d, s, , 1]
      v[s] <- if (sum(!is.na(vals)) >= 3) theil_sen_slope(vals) else NA_real_
    }
    v
  }, mc.cores = mc_cores)

  td_mat <- do.call(rbind, td_list)  # [draw, site]
  td_out <- t(td_mat)                 # [site, draw]  ← rf_importance_draws 期望格式
  rownames(td_out) <- sites
  colnames(td_out) <- NULL

  out_path <- v3_file("derived",
    paste0("trend_draws_", m, "_", run_label_ext), "rds")
  saveRDS(td_out, out_path)

  t_end <- Sys.time()
  message(sprintf("[05e]    saved %s (dim %d x %d, %.1fs)",
    basename(out_path), nrow(td_out), ncol(td_out),
    as.numeric(t_end - t_start, units = "secs")))
}

log_time("05e", "DONE")
