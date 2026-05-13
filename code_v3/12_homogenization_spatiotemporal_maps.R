#!/usr/bin/env Rscript
## 12_homogenization_spatiotemporal_maps.R  —  v3 均质化时空图
##
## 修复 v2 中 dist() 对经纬度用欧氏距离的 bug
## 使用测地线距离
## Mann-Kendall 正式检验
## v3 fix: 正确处理 4D psi（draws × species × sites × periods）
## v3 fix: 跨所有后验 draw 汇总（非仅第 1 个 draw）
## v3 fix: Sørensen 距离裁剪到 [0,1]

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(sf); library(Kendall)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_spatial.R"))
source(file.path(CODE_V3, "utils_mapping.R"))
P <- ensure_v3_dirs()

log_time("12", "Starting homogenization analysis")

is_pilot <- Sys.getenv("V3_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL

# ── 1. 加载数据 ──────────────────────────────────────────────────────
survey_stem <- paste0("survey_history", GRID_TAG, "_v3")
survey <- safe_read(v3_file("derived", survey_stem, "rds"))
if (is.null(survey)) {
  survey <- safe_read(v3_file("derived", "survey_history_v3", "rds"))
}
psi_obj <- safe_read(v3_file("derived", paste0("psi_samples_thinned_", run_label), "rds"))
grid_env_stem <- paste0("grid_environment", GRID_TAG, "_v3")
grid_env <- safe_read(v3_file("derived", grid_env_stem, "rds"))
if (is.null(grid_env)) {
  grid_env <- safe_read(v3_file("derived", "grid_environment_v3", "rds"))
}
if (is.null(grid_env)) {
  grid_env <- safe_read(file.path(DIRS$v2_derived, "grid_environment_dynamic_occupancy.rds"))
}

if (is.null(psi_obj)) stop("Psi samples not found. Run 05 first.")

psi_arr <- psi_obj$psi_samples_thinned
species <- psi_obj$species
sites   <- psi_obj$sites

# ── FIX #1: 正确推断 psi 维度和 period 数 ──────────────────────────
psi_dim <- dim(psi_arr)
if (length(psi_dim) == 4) {
  # 4D: draws × species × sites × periods
  n_draws  <- psi_dim[1]
  n_sp     <- psi_dim[2]
  n_sites  <- psi_dim[3]
  n_periods <- psi_dim[4]
  message(sprintf("[12] psi is 4D: %d draws × %d sp × %d sites × %d periods",
                  n_draws, n_sp, n_sites, n_periods))
} else if (length(psi_dim) == 3) {
  # 3D: draws × species × sites（单 period 或 period 被压平）
  n_draws  <- psi_dim[1]
  n_sp     <- psi_dim[2]
  n_sites  <- psi_dim[3]
  n_periods <- length(survey$periods)
  message(sprintf("[12] psi is 3D (%d × %d × %d), inferring n_periods=%d from survey",
                  n_draws, n_sp, n_sites, n_periods))
  if (n_periods == 1) {
    message("[12] Only 1 period available — homogenization analysis requires ≥2 periods. Exiting.")
    log_time("12", "SKIPPED (single period)")
    quit(save = "no", status = 0)
  }
} else {
  stop("[12] Unexpected psi dimensions: ", paste(psi_dim, collapse = "x"))
}

# ── 2. 计算群落相似性（Sørensen）─────────────────────────────────────
message("[12] Computing geodesic distance matrix ...")
geo_dist <- grid_distance_matrix(grid_env |> filter(grid_cell %in% sites))

# 取前几个后验抽取
n_draws_use <- min(40, n_draws)
draw_idx <- round(seq(1, n_draws, length.out = n_draws_use))

# 逐 draw × period 计算网格间 Sørensen 距离
sorensen_by_period <- list()
for (d in seq_len(n_draws_use)) {
  for (t in seq_len(n_periods)) {
    # FIX #1: 正确提取当前 draw + period 的 psi
    if (length(psi_dim) >= 4) {
      psi_dt <- psi_arr[draw_idx[d], , , t]  # species × sites
    } else {
      # 3D fallback: 无法区分 period，跳过
      next
    }
    # 二值化：psi > 0.5 视为存在
    pres <- (psi_dt > 0.5) * 1

    # Sørensen 距离矩阵（sites × sites）
    pres_t <- t(pres)  # sites × species
    a <- pres_t %*% t(pres_t)  # 共享物种数
    b <- rowSums(pres_t)       # 每个 site 的物种数
    denom <- outer(b, b, "+")
    # FIX #13: 裁剪到 [0,1]，防止浮点精度导致负值
    sorensen <- as.matrix(pmax(0, pmin(1, 1 - 2 * a / denom)))
    diag(sorensen) <- 0

    sorensen_by_period[[paste0("d", d, "_P", t)]] <- sorensen
  }
  if (d %% 10 == 0) message(sprintf("  [12] ... draw %d/%d", d, n_draws_use))
}

# ── 3. Mann-Kendall 趋势检验 ──────────────────────────────────────────
# FIX #2: 跨所有 draw 取平均（非仅第 1 个 draw）

mean_sorensen <- numeric(n_periods)
sorensen_sd   <- numeric(n_periods)
for (t in seq_len(n_periods)) {
  # 收集该 period 所有 draw 的上三角均值
  vals <- numeric(n_draws_use)
  for (d in seq_len(n_draws_use)) {
    key <- paste0("d", d, "_P", t)
    if (key %in% names(sorensen_by_period)) {
      s_mat <- sorensen_by_period[[key]]
      vals[d] <- mean(s_mat[upper.tri(s_mat)], na.rm = TRUE)
    } else {
      vals[d] <- NA_real_
    }
  }
  mean_sorensen[t] <- mean(vals, na.rm = TRUE)
  sorensen_sd[t]   <- sd(vals, na.rm = TRUE)
}

mk_test <- MannKendall(mean_sorensen)
message(sprintf("[12] Mann-Kendall test for mean Sørensen distance: tau = %.3f, p = %.4f",
                mk_test$tau[1], mk_test$sl[1]))

# ── 4. 均质化趋势图（带后验不确定性）─────────────────────────────────
mk_df <- tibble(
  period = paste0("P", seq_len(n_periods)),
  mean_sorensen = mean_sorensen,
  sd_sorensen = sorensen_sd
)

p_homog <- ggplot(mk_df, aes(x = period, y = mean_sorensen)) +
  geom_ribbon(aes(ymin = mean_sorensen - sd_sorensen,
                  ymax = mean_sorensen + sd_sorensen),
              alpha = 0.15, fill = NATURE_ACCENT) +
  geom_line(colour = NATURE_ACCENT, linewidth = 0.5, group = 1) +
  geom_point(colour = NATURE_ACCENT, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey50",
               linewidth = 0.3, fill = "grey80", alpha = 0.3) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.3,
           label = sprintf("Mann-Kendall τ = %.2f, p = %.3f",
                           mk_test$tau[1], mk_test$sl[1]),
           size = 2.5, colour = "grey30") +
  labs(x = "Period", y = "Mean Sørensen distance",
       title = "Biotic homogenization trend") +
  theme_nature_pub()

save_nature(p_homog, paste0("fig_homogenization_trend_", run_label),
             width_mm = 89, height_mm = 60)

# ── 5. 保存结果 ──────────────────────────────────────────────────────
write_csv(mk_df, v3_file("results", paste0("table_homogenization_trend_", run_label)))

log_time("12", "DONE")
