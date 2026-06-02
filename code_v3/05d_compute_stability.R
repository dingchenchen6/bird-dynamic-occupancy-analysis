#!/usr/bin/env Rscript
## 05d_compute_stability.R  —  群落稳定性、韧性与抗性分析
##
## 输入: metric_arrays_*.rds (由 05_postprocess_diversity_extended.R 生成)
##       grid_environment_v3.rds (环境/扰动数据)
## 输出: table_stability_summary_*.csv
##       table_resistance_summary_*.csv

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
P <- ensure_v3_dirs()

run_label_extended <- paste0(RUN_LABEL, "_extended")
log_time("05d", "Starting stability analysis")

# ── 1. 加载数据 ──────────────────────────────────────────────────────
metric_data <- safe_read(v3_file("derived", paste0("metric_arrays_", run_label_extended), "rds"))
if (is.null(metric_data)) stop("Metric arrays not found. Run 05_postprocess_diversity_extended.R first.")

metric_arr    <- metric_data$metric_arr      # draws × sites × periods × metrics
temporal_arr  <- metric_data$temporal_arr    # draws × sites × temporal_metrics
species       <- metric_data$species
sites         <- metric_data$sites
n_periods     <- metric_data$n_periods
metrics       <- metric_data$metrics
n_draws       <- dim(metric_arr)[1]
n_sites       <- dim(metric_arr)[2]

grid_env <- safe_read(v3_file("derived", paste0("grid_environment", GRID_TAG, "_v3"), "rds"))
if (is.null(grid_env)) {
  grid_env <- safe_read(v3_file("derived", "grid_environment_v3", "rds"))
}

# ── 2. 时间稳定性指标 ────────────────────────────────────────────────
# CV^{-1} = mean / sd (越高越稳定)
# 对 5 个 period 的每个 metric，计算每个 grid 的 CV^{-1}

stability_metrics <- c("temporal_stability_cv_inv", "temporal_stability_logcv")
stability_arr <- array(NA_real_,
                       dim = c(n_draws, n_sites, length(metrics), length(stability_metrics)),
                       dimnames = list(draw = NULL, site = sites,
                                       metric = metrics, stability = stability_metrics))

message(sprintf("[05d] Computing temporal stability for %d draws × %d sites × %d metrics",
                n_draws, n_sites, length(metrics)))

for (d in seq_len(n_draws)) {
  for (s in seq_len(n_sites)) {
    for (m in metrics) {
      vals <- metric_arr[d, s, , m]
      if (sum(!is.na(vals)) >= 3) {
        mmean <- mean(vals, na.rm = TRUE)
        msd   <- sd(vals, na.rm = TRUE)
        if (msd > 0) {
          stability_arr[d, s, m, "temporal_stability_cv_inv"] <- mmean / msd
          stability_arr[d, s, m, "temporal_stability_logcv"]  <- log(msd / mmean)
        }
      }
    }
  }
  if (d %% 50 == 0) message(sprintf("  [05d] ... draw %d/%d", d, n_draws))
}

# ── 3. 组成稳定性：轨迹距离 ─────────────────────────────────────────
# 对每个 draw × grid，计算物种组成从 P1 到 P5 的 Euclidean 距离
# 以及相邻 period 间的平均 Bray-Curtis 距离

composition_stability_arr <- array(NA_real_, dim = c(n_draws, n_sites, 2),
  dimnames = list(draw = NULL, site = sites,
                  metric = c("trajectory_distance", "mean_bray_curtis_change")))

message("[05d] Computing composition stability ...")

# 需要 psi_samples 才能计算组成稳定性
# 如果 psi_samples 不在内存中，尝试读取 thinned psi
psi_path <- v3_file("derived", paste0("psi_samples_thinned_", run_label_extended), "rds")
psi_data <- safe_read(psi_path)
if (is.null(psi_data)) {
  warning("[05d] psi_samples_thinned not found. Skipping composition stability.")
} else {
  psi_samples <- psi_data$psi_samples_thinned
  psi_dim <- dim(psi_samples)
  draw_idx <- psi_data$draw_idx

  for (d in seq_len(n_draws)) {
    for (s in seq_len(n_sites)) {
      if (length(psi_dim) >= 4) {
        psi_site <- psi_samples[d, , s, ]  # species × periods
      } else {
        next
      }
      if (!is.matrix(psi_site)) next
      if (ncol(psi_site) < 2) next

      # 轨迹距离：P1 到 P5 的 Euclidean 距离
      composition_stability_arr[d, s, "trajectory_distance"] <-
        sqrt(sum((psi_site[, 1] - psi_site[, ncol(psi_site)])^2))

      # 平均 Bray-Curtis 变化
      bc_changes <- numeric(ncol(psi_site) - 1)
      for (t in seq_len(ncol(psi_site) - 1)) {
        p1 <- psi_site[, t]; p2 <- psi_site[, t + 1]
        bc_changes[t] <- sum(abs(p1 - p2)) / sum(p1 + p2)
      }
      composition_stability_arr[d, s, "mean_bray_curtis_change"] <- mean(bc_changes)
    }
    if (d %% 50 == 0) message(sprintf("  [05d composition] ... draw %d/%d", d, n_draws))
  }
}

# ── 4. 韧性 (Resilience)：AR1 自回归系数 ─────────────────────────────
# 对多样性异常的 AR1：anomaly_t = phi * anomaly_{t-1} + epsilon
# phi 越接近 0 = 回归越快 = 韧性越强

resilience_arr <- array(NA_real_, dim = c(n_draws, n_sites, length(metrics)),
  dimnames = list(draw = NULL, site = sites, metric = metrics))

message("[05d] Computing resilience (AR1) ...")

for (d in seq_len(n_draws)) {
  for (s in seq_len(n_sites)) {
    for (m in metrics) {
      vals <- metric_arr[d, s, , m]
      if (sum(!is.na(vals)) >= 4) {
        # 计算异常值（去趋势）
        trend <- coef(lm(vals ~ seq_len(n_periods)))[2]
        anomaly <- vals - (mean(vals, na.rm = TRUE) + trend * (seq_len(n_periods) - (n_periods + 1) / 2))
        # AR1: anomaly[t] ~ anomaly[t-1]
        if (var(anomaly[-1], na.rm = TRUE) > 0) {
          ar1_fit <- tryCatch(
            lm(anomaly[-1] ~ anomaly[-n_periods] - 1),
            error = function(e) NULL
          )
          if (!is.null(ar1_fit)) {
            resilience_arr[d, s, m] <- coef(ar1_fit)[1]
          }
        }
      }
    }
  }
  if (d %% 50 == 0) message(sprintf("  [05d resilience] ... draw %d/%d", d, n_draws))
}

# ── 5. 抗性 (Resistance)：扰动期变化幅度 ────────────────────────────
# 使用 HFI 变化 (delta_hfi) 和气温变化 (delta_t_mean) 作为扰动代理
# Resistance = |Delta diversity during high perturbation| / |Delta during low perturbation|

resistance_arr <- array(NA_real_, dim = c(n_draws, n_sites, length(metrics), 2),
  dimnames = list(draw = NULL, site = sites, metric = metrics,
                  perturbation = c("hfi_change", "climate_change")))

if (!is.null(grid_env) && "delta_hfi" %in% names(grid_env)) {
  message("[05d] Computing resistance ...")
  # 简化：使用 grid_env 中的 delta 变量（代表 P1→P5 的总变化）
  # 高扰动 = delta 绝对值上三分位数
  env_sub <- grid_env |> filter(grid_cell %in% sites)
  env_sub$delta_hfi_rank <- rank(abs(env_sub$delta_hfi))
  high_perturb_hfi <- env_sub$delta_hfi_rank > quantile(env_sub$delta_hfi_rank, 0.67, na.rm = TRUE)

  for (d in seq_len(n_draws)) {
    for (m in metrics) {
      vals_start <- metric_arr[d, , 1, m]
      vals_end   <- metric_arr[d, , n_periods, m]
      delta_div <- abs(vals_end - vals_start)

      if (sum(high_perturb_hfi, na.rm = TRUE) > 0) {
        resistance_arr[d, , m, "hfi_change"] <- delta_div
      }
    }
  }
} else {
  message("[05d] grid_env delta variables not found. Skipping resistance.")
}

# ── 6. 汇总后验并输出 ───────────────────────────────────────────────

summarise_stability <- function(arr, metric_name, extra_cols = NULL) {
  n_draws_ <- dim(arr)[1]
  n_sites_ <- dim(arr)[2]
  out <- tibble()
  for (s in seq_len(n_sites_)) {
    vals <- arr[, s]
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) next
    out <- bind_rows(out, tibble(
      grid_cell = sites[s],
      metric = metric_name,
      mean = mean(vals),
      sd = sd(vals),
      q025 = quantile(vals, 0.025),
      median = median(vals),
      q975 = quantile(vals, 0.975),
      !!!extra_cols
    ))
  }
  out
}

# 时间稳定性汇总
stability_summary <- list()
for (m in metrics) {
  for (st in stability_metrics) {
    smry <- summarise_stability(stability_arr[, , m, st], m, list(stability_metric = st))
    stability_summary[[length(stability_summary) + 1]] <- smry
  }
}
# 组成稳定性
for (st in c("trajectory_distance", "mean_bray_curtis_change")) {
  smry <- summarise_stability(composition_stability_arr[, , st], "composition", list(stability_metric = st))
  stability_summary[[length(stability_summary) + 1]] <- smry
}
# 韧性
for (m in metrics) {
  smry <- summarise_stability(resilience_arr[, , m], m, list(stability_metric = "ar1_resilience"))
  stability_summary[[length(stability_summary) + 1]] <- smry
}

stability_df <- bind_rows(stability_summary)
write_csv(stability_df, v3_file("results", paste0("table_stability_summary_", run_label_extended)))
message(sprintf("[05d] Stability summary: %d rows", nrow(stability_df)))

# 抗性汇总
resistance_summary <- list()
for (m in metrics) {
  for (pert in c("hfi_change", "climate_change")) {
    smry <- summarise_stability(resistance_arr[, , m, pert], m, list(perturbation = pert))
    resistance_summary[[length(resistance_summary) + 1]] <- smry
  }
}
resistance_df <- bind_rows(resistance_summary)
if (nrow(resistance_df) > 0) {
  write_csv(resistance_df, v3_file("results", paste0("table_resistance_summary_", run_label_extended)))
  message(sprintf("[05d] Resistance summary: %d rows", nrow(resistance_df)))
}

log_time("05d", "DONE: stability analysis")
message("[05d] All stability metrics computed and saved.")
