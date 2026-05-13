#!/usr/bin/env Rscript
## 05_postprocess_diversity.R  —  v3 后处理与多样性计算
##
## 关键修复：
##   1. log10 Bug: if_else(.x > 0, log10(.x), NA_real_)
##   2. brms 加 gp(centroid_lon, centroid_lat) 空间高斯过程
##   3. 使用 trait_extended_v3.rds（含 habitat_breadth + diet_specialization）
##   4. adapt_delta=0.99, max_treedepth=15
##   5. centroid_lon 加入驱动变量
##   6. v3 新增 codyn 风格时间动态指标（synchrony / variance_ratio / turnover）
##   7. 修复 n_periods 推断：从 psi_samples 维度推导（非 length(fit$y)）

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(brms); library(cmdstanr); library(abind)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_diversity.R"))
P <- ensure_v3_dirs()

log_time("05", "Starting postprocess diversity")

# ── 1. 加载模型 ──────────────────────────────────────────────────────
is_pilot <- Sys.getenv("V3_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL

# 兼容 tMsPGOcc 和 stMsPGOcc 输出文件名
fit <- safe_read(v3_file("derived", paste0("stMsPGOcc_fit_", run_label), "rds"))
if (is.null(fit)) {
  fit <- safe_read(v3_file("derived", paste0("tMsPGOcc_fit_", run_label), "rds"))
}
if (is.null(fit)) stop("Model fit not found (tried stMsPGOcc and tMsPGOcc). Run 04 first.")

grid_env_stem <- paste0("grid_environment", GRID_TAG, "_v3")
grid_env <- safe_read(v3_file("derived", grid_env_stem, "rds"))
if (is.null(grid_env)) {
  grid_env <- safe_read(v3_file("derived", "grid_environment_v3", "rds"))
}
if (is.null(grid_env)) {
  grid_env <- safe_read(file.path(DIRS$v2_derived, "grid_environment_dynamic_occupancy.rds"))
}

# ── 2. 提取 psi 后验 ──────────────────────────────────────────────────
psi_samples <- fit$psi.samples  # draws × species × sites
n_draws <- min(dim(psi_samples)[1], PSI_MAX_DRAWS)
draw_idx <- round(seq(1, dim(psi_samples)[1], length.out = n_draws))

n_sp <- dim(psi_samples)[2]
n_sites <- dim(psi_samples)[3]
species <- rownames(fit$y) %||% paste0("sp", seq_len(n_sp))
sites <- colnames(fit$y) %||% paste0("site", seq_len(n_sites))

# ── 3. 加载扩展性状（含 habitat_breadth + diet_specialization）───────
trait_ext <- safe_read(v3_file("derived", "trait_extended_v3", "rds"))
if (is.null(trait_ext)) {
  warning("[05] trait_extended_v3 not found. Falling back to basic traits.")
  trait_ext <- read_csv_safe(TRAIT_IMPUTED_PATH) |>
    mutate(diet_specialization = NA_real_, habitat_breadth = NA_real_)
}

# 使用修复后的 prepare_trait_matrix（v2 log10 bug 已修复）
trait_aligned <- tibble(species = as.character(species)) |>
  left_join(trait_ext %>% mutate(species = as.character(species)), by = "species")

# 合并基础 + 扩展性状
all_trait_vars <- intersect(
  TRAIT_VARS_ALL,
  names(trait_aligned)
)
trait_mat <- prepare_trait_matrix(trait_aligned, all_trait_vars)

# 距离矩阵
dist_mat <- as.matrix(stats::dist(trait_mat))

# ── 4. 加载系统发育 ──────────────────────────────────────────────────
phylo <- NULL
if (requireNamespace("ape", quietly = TRUE)) {
  phylo_path <- v3_file("derived", "phylogeny_matched_v3", "rds")
  if (!file.exists(phylo_path)) {
    phylo_path <- file.path(DIRS$v2_derived, "phylogeny_matched.rds")
  }
  phylo <- safe_read(phylo_path)
}

# ── 5. 计算多样性指标（传播后验不确定性）───────────────────────────────
# 修复 n_periods 推断：从 psi_samples 维度推导
# stMsPGOcc: psi_samples 为 4D (draws × species × sites × periods)
# msPGOcc:   psi_samples 为 3D (draws × species × sites)，仅 1 个 period
psi_dim <- dim(psi_samples)
n_periods <- if (length(psi_dim) >= 4) psi_dim[4] else 1L

metrics <- c("corrected_richness", "shannon", "inv_simpson",
              "trait_volume", "trait_dispersion", "rao_q",
              "feve", "fdiv",
              "pd_prob", "mpd_prob")

# 创建存储数组
metric_arr <- array(NA_real_,
                     dim = c(n_draws, n_sites, n_periods, length(metrics)),
                     dimnames = list(
                       draw = NULL,
                       site = sites,
                       period = NULL,
                       metric = metrics
                     ))

message(sprintf("[05] Computing diversity for %d draws × %d sites × %d periods",
                n_draws, n_sites, n_periods))

for (d in seq_len(n_draws)) {
  for (t in seq_len(n_periods)) {
    # 从 psi_samples 正确提取当前 draw + period 的 psi 矩阵
    # stMsPGOcc 4D: psi_samples[draw, species, sites, period]
    # msPGOcc   3D: psi_samples[draw, species, sites]（period=1）
    if (length(psi_dim) >= 4) {
      psi_d <- psi_samples[draw_idx[d], , , t]  # species × sites
    } else {
      psi_d <- psi_samples[draw_idx[d], , ]      # species × sites
    }
    if (is.matrix(psi_d)) {
      for (s in seq_len(n_sites)) {
        psi_v <- psi_d[, s]

        # 分类学多样性
        tax <- div_taxonomic(psi_v)
        metric_arr[d, s, t, "corrected_richness"] <- tax$richness
        metric_arr[d, s, t, "shannon"]             <- tax$shannon
        metric_arr[d, s, t, "inv_simpson"]          <- tax$inv_simpson

        # 功能多样性
        fn <- div_functional(psi_v, trait_mat)
        metric_arr[d, s, t, "trait_volume"]     <- fn$trait_volume
        metric_arr[d, s, t, "trait_dispersion"] <- fn$trait_dispersion
        metric_arr[d, s, t, "rao_q"]            <- fn$rao_q

        # 功能均匀度 FEve（Villéger et al. 2008）
        metric_arr[d, s, t, "feve"] <- div_feve(psi_v, trait_mat)

        # 功能分散度 FDiv（Mason et al. 2003; Villéger et al. 2008）
        metric_arr[d, s, t, "fdiv"] <- div_fdiv(psi_v, trait_mat)

        # 系统发育多样性
        if (!is.null(phylo)) {
          pd <- div_phylogenetic(psi_v, phylo, tip_order = species)
          metric_arr[d, s, t, "pd_prob"]  <- pd$pd_prob
          metric_arr[d, s, t, "mpd_prob"] <- pd$mpd_prob
        }
      }
    }
  }
  if (d %% 50 == 0) message(sprintf("  ... draw %d/%d", d, n_draws))
}

# ── 5b. 计算群落时间动态指标（codyn 风格）──────────────────────────
# 对每个网格 × 每个后验 draw，构建 species × periods 群落矩阵
# 计算：同步性 (synchrony)、方差比 (variance_ratio)、时间周转 (turnover)
# 这些是站点级指标（跨所有 period 计算），与逐 period 的多样性指标不同

temporal_metrics <- c("synchrony", "variance_ratio",
                       "turnover_total", "turnover_gain", "turnover_loss")

temporal_arr <- array(NA_real_,
                       dim = c(n_draws, n_sites, length(temporal_metrics)),
                       dimnames = list(draw = NULL, site = sites, metric = temporal_metrics))

message(sprintf("[05] Computing temporal dynamics for %d draws × %d sites",
                n_draws, n_sites))

for (d in seq_len(n_draws)) {
  for (s in seq_len(n_sites)) {
    # 构建 species × periods 群落矩阵
    if (length(psi_dim) >= 4) {
      comm_mat <- psi_samples[draw_idx[d], , s, ]  # species × periods
    } else {
      # 静态模型只有 1 个 period → 无法计算时间动态
      next
    }

    # 确保是矩阵（species × periods）
    if (!is.matrix(comm_mat)) {
      comm_mat <- matrix(comm_mat, nrow = n_sp, ncol = 1)
    }

    if (ncol(comm_mat) < 2) next  # 需要至少 2 个 period

    # 同步性
    temporal_arr[d, s, "synchrony"] <- div_synchrony(comm_mat)

    # 方差比
    temporal_arr[d, s, "variance_ratio"] <- div_variance_ratio(comm_mat)

    # 时间周转
    to <- div_temporal_turnover(comm_mat, threshold = 0.1)
    temporal_arr[d, s, "turnover_total"] <- to$turnover_total
    temporal_arr[d, s, "turnover_gain"]  <- to$turnover_gain
    temporal_arr[d, s, "turnover_loss"]  <- to$turnover_loss
  }
  if (d %% 50 == 0) message(sprintf("  [temporal] ... draw %d/%d", d, n_draws))
}

message(sprintf("[05] Temporal dynamics computed. Non-NA rates: %s",
  paste(temporal_metrics, "=",
        round(colSums(!is.na(temporal_arr[, , ])) / (n_draws * n_sites) * 100, 1),
        collapse = ", ")))

# ── 6. 汇总后验 ──────────────────────────────────────────────────────
# 网格 × period × metric 的后验均值和 95% CRI
metric_summary <- list()
for (m in metrics) {
  arr_m <- metric_arr[, , , m]
  for (t in seq_len(n_periods)) {
    smry <- summarise_post(arr_m[, , t])
    smry$metric <- m
    smry$period <- paste0("P", t)
    smry$grid_cell <- sites
    metric_summary[[length(metric_summary) + 1]] <- smry
  }
}
metric_summary_df <- bind_rows(metric_summary)
write_csv(metric_summary_df, v3_file("results", paste0("table_diversity_summary_", run_label)))

# 同时输出宽格式（兼容06_figures_publication.R）
metric_wide <- metric_summary_df |>
  pivot_wider(
    id_cols = c(grid_cell, period),
    names_from = metric,
    values_from = c(mean, sd, q025, q975),
    names_glue = "{metric}_{.value}"
  )
write_csv(metric_wide, v3_file("results", paste0("table_diversity_wide_", run_label)))

# 输出06需要的长格式（value_mean, value_l95, value_u95）
metric_cri_long <- metric_summary_df |>
  select(grid_cell, period, metric, value_mean = mean, value_l95 = q025, value_u95 = q975) |>
  mutate(block_label = period)
write_csv(metric_cri_long, v3_file("results", paste0("table_community_metrics_with_cri_", run_label)))

# ── 6b. 汇总时间动态后验 ────────────────────────────────────────────
# 网格 × temporal_metric 的后验均值和 95% CRI（无 period 维度）
temporal_summary <- list()
for (m in temporal_metrics) {
  smry <- summarise_post(temporal_arr[, , m])
  smry$metric <- m
  smry$grid_cell <- sites
  temporal_summary[[length(temporal_summary) + 1]] <- smry
}
temporal_summary_df <- bind_rows(temporal_summary)
write_csv(temporal_summary_df, v3_file("results", paste0("table_temporal_dynamics_summary_", run_label)))
message(sprintf("[05] Temporal dynamics summary: %d rows × %d cols, %d metrics",
               nrow(temporal_summary_df), ncol(temporal_summary_df), length(temporal_metrics)))

# ── 5c. 概率加权 Baselga 时间分解 ──────────────────────────────────────
# Q2 核心：β 多样性中周转 vs 嵌套的比例
# 对每个网格 × 每个后验 draw，计算相邻时期对的 Baselga 分解
# 输出：beta_sor, beta_sim (turnover), beta_sne (nestedness), prop_turnover

baselga_metrics <- c("beta_sor", "beta_sim", "beta_sne", "prop_turnover")
n_period_pairs <- n_periods - 1

if (n_period_pairs >= 1) {
  baselga_arr <- array(NA_real_,
                       dim = c(n_draws, n_sites, n_period_pairs, length(baselga_metrics)),
                       dimnames = list(draw = NULL, site = sites,
                                       period_pair = NULL, metric = baselga_metrics))

  message(sprintf("[05] Computing Baselga decomposition for %d draws × %d sites × %d period-pairs",
                  n_draws, n_sites, n_period_pairs))

  for (d in seq_len(n_draws)) {
    for (s in seq_len(n_sites)) {
      # 构建 species × periods 群落矩阵
      if (length(psi_dim) >= 4) {
        comm_mat <- psi_samples[draw_idx[d], , s, ]  # species × periods
      } else {
        next
      }
      if (!is.matrix(comm_mat)) comm_mat <- matrix(comm_mat, nrow = n_sp, ncol = n_periods)
      if (ncol(comm_mat) < 2) next

      for (pp in seq_len(n_period_pairs)) {
        bg <- div_baselga(comm_mat[, pp], comm_mat[, pp + 1])
        baselga_arr[d, s, pp, "beta_sor"]       <- bg$beta_sor
        baselga_arr[d, s, pp, "beta_sim"]       <- bg$beta_sim
        baselga_arr[d, s, pp, "beta_sne"]       <- bg$beta_sne
        baselga_arr[d, s, pp, "prop_turnover"]  <- bg$prop_turnover
      }
    }
    if (d %% 50 == 0) message(sprintf("  [baselga] ... draw %d/%d", d, n_draws))
  }

  # 汇总 Baselga 后验
  baselga_summary <- list()
  for (pp in seq_len(n_period_pairs)) {
    for (m in baselga_metrics) {
      smry <- summarise_post(baselga_arr[, , pp, m])
      smry$metric <- m
      smry$period_pair <- paste0("P", pp, "_P", pp + 1)
      smry$grid_cell <- sites
      baselga_summary[[length(baselga_summary) + 1]] <- smry
    }
  }
  baselga_summary_df <- bind_rows(baselga_summary)
  write_csv(baselga_summary_df, v3_file("results", paste0("table_baselga_summary_", run_label)))

  # 全局汇总：周转占比的跨网格均值
  global_baselga <- baselga_summary_df |>
    filter(metric == "prop_turnover") |>
    group_by(period_pair) |>
    summarise(
      prop_turnover_mean = mean(mean, na.rm = TRUE),
      prop_turnover_q025 = mean(q025, na.rm = TRUE),
      prop_turnover_q975 = mean(q975, na.rm = TRUE),
      n_grids = sum(!is.na(mean)),
      .groups = "drop"
    )
  write_csv(global_baselga, v3_file("results", paste0("table_baselga_global_", run_label)))

  message(sprintf("[05] Baselga: prop_turnover (global mean) = %s",
    paste(global_baselga$period_pair, "=",
          round(global_baselga$prop_turnover_mean, 3), collapse = ", ")))
} else {
  message("[05] Skipping Baselga decomposition (only 1 period)")
}

# ── 7. 计算趋势（FIX #6: OLS + Theil-Sen 双重估计）──────────────────
trend_methods <- c("ols", "theil_sen")
trend_draws <- array(NA_real_,
                       dim = c(n_draws, n_sites, length(metrics), length(trend_methods)),
                       dimnames = list(draw = NULL, site = sites, metric = metrics, method = trend_methods))

for (d in seq_len(n_draws)) {
  for (s in seq_len(n_sites)) {
    for (m in metrics) {
      vals <- metric_arr[d, s, , m]
      if (sum(!is.na(vals)) >= 3) {
        # OLS 斜率（保留原有方法）
        trend_draws[d, s, m, "ols"] <- coef(lm(vals ~ seq_len(n_periods)))[2]
        # Theil-Sen 斜率（对端点异常值不敏感）
        trend_draws[d, s, m, "theil_sen"] <- theil_sen_slope(vals)
      }
    }
  }
}

# 趋势后验汇总（两种方法分别输出）
trend_summary <- list()
for (method in trend_methods) {
  for (m in metrics) {
    smry <- summarise_post(trend_draws[, , m, method])
    smry$metric <- m
    smry$method <- method
    smry$grid_cell <- sites
    trend_summary[[length(trend_summary) + 1]] <- smry
  }
}
trend_summary_df <- bind_rows(trend_summary)
write_csv(trend_summary_df, v3_file("results", paste0("table_trend_summary_", run_label)))

# Mann-Kendall 检验（基于后验均值的非参数趋势检验，补充 OLS 的 p 值）
mk_results <- list()
for (s in seq_len(n_sites)) {
  for (m in metrics) {
    # 取后验均值时间序列
    vals_mean <- apply(metric_arr[, s, , m], 2, mean, na.rm = TRUE)
    if (sum(!is.na(vals_mean)) >= 4) {
      mk <- mann_kendall_test(vals_mean)
      mk_results[[length(mk_results) + 1]] <- tibble(
        grid_cell = sites[s], metric = m,
        mk_S = mk$S, mk_tau = mk$tau, mk_p = mk$p_value
      )
    }
  }
}
mk_df <- bind_rows(mk_results)
write_csv(mk_df, v3_file("results", paste0("table_mann_kendall_", run_label)))
message(sprintf("[05] Mann-Kendall test: %d grid-metric combinations tested", nrow(mk_df)))

# ── 7b. Naive vs Corrected 对比（FIX #16）─────────────────────────────
# Q5 核心问题：占域校正改变了哪些结论？
# 计算：(1) 趋势方向翻转比例 (2) 趋势量级偏差 (3) 配对比较

# 朴素丰富度 = 每个 period 的原始物种计数（不考虑检测概率）
# 从原始 y 矩阵计算
if (exists("fit") && !is.null(fit$y)) {
  message("[05] Computing naive vs corrected richness comparison")

  # y 是 species × site_period 矩阵或 4D array
  y_data <- fit$y

  # 朴素丰富度趋势（基于原始出现/缺失）
  if (is.list(y_data) && length(y_data) == n_periods) {
    # y 为 period 列表
    naive_richness <- sapply(y_data, function(mat) colSums(mat, na.rm = TRUE))
  } else if (is.array(y_data) && length(dim(y_data)) == 4) {
    # y 为 4D array (species × sites × primary × secondary)
    # FIX: 先跨 secondary 合并（任一 rep 检测到 = 存在），再跨 species 计数
    det_merged <- apply(y_data > 0, c(1, 2, 3), any, na.rm = TRUE)
    naive_richness <- apply(det_merged, c(2, 3), sum, na.rm = TRUE)
  } else {
    naive_richness <- NULL
  }

  if (!is.null(naive_richness) && ncol(naive_richness) == n_periods) {
    # 校正丰富度后验均值
    corrected_richness_mean <- apply(metric_arr[, , , "corrected_richness"], c(2, 3), mean, na.rm = TRUE)

    # 计算趋势
    naive_trends <- apply(naive_richness, 1, function(x) {
      if (sum(!is.na(x)) >= 3) coef(lm(x ~ seq_len(n_periods)))[2] else NA_real_
    })
    corrected_trends <- apply(corrected_richness_mean, 1, function(x) {
      if (sum(!is.na(x)) >= 3) coef(lm(x ~ seq_len(n_periods)))[2] else NA_real_
    })

    # 趋势方向翻转
    naive_sign <- sign(naive_trends)
    corrected_sign <- sign(corrected_trends)
    flip_idx <- which(naive_sign != corrected_sign & !is.na(naive_sign) & !is.na(corrected_sign) & naive_sign != 0 & corrected_sign != 0)
    flip_rate <- length(flip_idx) / sum(!is.na(naive_sign) & !is.na(corrected_sign) & naive_sign != 0 & corrected_sign != 0)

    # 趋势量级偏差
    trend_diff <- corrected_trends - naive_trends
    trend_diff_abs <- abs(trend_diff[!is.na(trend_diff)])

    comparison_df <- tibble(
      grid_cell = sites,
      naive_trend = naive_trends,
      corrected_trend = corrected_trends,
      trend_diff = trend_diff,
      direction_flipped = naive_sign != corrected_sign
    )

    write_csv(comparison_df, v3_file("results", paste0("table_naive_vs_corrected_", run_label)))

    message(sprintf("[05] Naive vs Corrected: flip_rate = %.1f%%, median |trend_diff| = %.3f",
                   flip_rate * 100, median(trend_diff_abs, na.rm = TRUE)))
  } else {
    message("[05] Cannot compute naive richness (y format not recognized)")
  }
} else {
  message("[05] Skipping naive vs corrected comparison (fit$y not available)")
}

# ── 8. brms 驱动回归（加空间 GP）────────────────────────────────────────
# 对每种指标的趋势拟合 brms
driver_vars <- c("bio4", "bio7", "bio11", "bio13",
                  "elev_mean", "elev_sd", "texture_shannon",
                  "habitat_diversity_shannon",
                  "hfi_mean", "landcover_built", "landcover_cropland",
                  "centroid_lon", "centroid_lat")

for (m in c("corrected_richness", "shannon", "trait_volume", "pd_prob")) {
  message(sprintf("[05] brms driver regression for: %s", m))

  trend_m <- trend_summary_df |>
    filter(metric == m) |>
    inner_join(grid_env, by = "grid_cell") |>
    drop_na()

  # 标准化驱动变量
  for (v in driver_vars) {
    if (v %in% names(trend_m)) {
      trend_m[[paste0("z_", v)]] <- scale(trend_m[[v]])[, 1]
    }
  }

  # brms 公式：加 gp(centroid_lon, centroid_lat) 空间高斯过程
  z_vars <- paste0("z_", driver_vars)
  z_vars <- intersect(z_vars, names(trend_m))
  # 移除 z_centroid_lon 和 z_centroid_lat（已在 gp 中）
  z_vars_gp <- setdiff(z_vars, c("z_centroid_lon", "z_centroid_lat"))

  brms_formula <- as.formula(
    paste0("mean ~ ", paste(z_vars_gp, collapse = " + "),
           " + gp(centroid_lon, centroid_lat, k = 10)")
  )

  brms_prior <- c(
    prior(normal(0, 1), class = "b"),
    prior(normal(0, 1), class = "Intercept")
  )

  brms_fit <- tryCatch({
    brm(
      formula    = brms_formula,
      data       = trend_m,
      family     = gaussian(),
      prior      = brms_prior,
      iter       = BRMS_ITER,
      warmup     = BRMS_WARMUP,
      chains     = BRMS_CHAINS,
      cores      = max(1, parallel::detectCores() - 1),
      control    = list(adapt_delta = BRMS_ADAPT_DELTA,
                         max_treedepth = BRMS_MAX_TREED),
      seed       = BRMS_SEED,
      backend    = "cmdstanr",
      silent     = 2
    )
  }, error = function(e) {
    warning(sprintf("[05] brms failed for %s: %s", m, e$message))
    NULL
  })

  if (!is.null(brms_fit)) {
    checkpoint_save(brms_fit, paste0("brms_driver_trend_", m, "_", run_label),
                     subdir = "derived")

    # DHARMa 检查
    dharma_res <- check_dharma_gate(
      brms_fit,
      save_dir = DIRS$figures,
      stem = paste0("dharma_driver_trend_", m, "_", run_label)
    )

    # LOO
    loo_res <- tryCatch(lo(brms_fit), error = function(e) NULL)
    if (!is.null(loo_res)) {
      checkpoint_save(loo_res, paste0("brms_loo_trend_", m, "_", run_label),
                       subdir = "derived")
    }
  }
}

# ── 9. 保存 psi 抽取 ──────────────────────────────────────────────────
# FIX #3: 添加维度断言，确保 psi_samples_thinned 保存时维度正确
if (length(psi_dim) >= 4) {
  psi_thinned <- psi_samples[draw_idx, , , , drop = FALSE]
} else {
  psi_thinned <- psi_samples[draw_idx, , , drop = FALSE]
}
psi_thinned_dim <- dim(psi_thinned)
expected_dim <- if (length(psi_dim) >= 4) {
  c(length(draw_idx), n_sp, n_sites, n_periods)
} else {
  c(length(draw_idx), n_sp, n_sites)
}
if (length(psi_thinned_dim) != length(expected_dim) ||
    any(psi_thinned_dim != expected_dim)) {
  warning(sprintf(
    "[05] psi_thinned dimension mismatch: got [%s], expected [%s]. Saving anyway.",
    paste(psi_thinned_dim, collapse = "x"),
    paste(expected_dim, collapse = "x")))
}
saveRDS(
  list(psi_samples_thinned = psi_thinned,
       species = species, sites = sites, draw_idx = draw_idx,
       psi_dim = psi_thinned_dim, n_periods = n_periods),
  v3_file("derived", paste0("psi_samples_thinned_", run_label), "rds")
)

# ═══════════════════════════════════════════════════════════════════════
# ── 10. 物种水平趋势系统分析（Q2/Q4 核心故事）────────────────────────
# ═══════════════════════════════════════════════════════════════════════
# 哪些物种在扩张？哪些在收缩？与气候变化、土地利用变化、性状有何关系？

message("[05] Starting species-level trend analysis (Section 10)")

# ── 10a. 计算每个物种的时间趋势 ──────────────────────────────────────
# 对每个后验 draw，计算每物种的 psi 时间趋势（跨所有 site 取加权均值）
# 趋势方法：Theil-Sen（稳健）+ OLS
n_sp_use <- length(species)

sp_trend_arr <- array(NA_real_,
                       dim = c(n_draws, n_sp_use, 2),
                       dimnames = list(draw = NULL, species = species,
                                       method = c("theil_sen", "ols")))

# 同时记录每物种每 period 的全空间均值 psi（用于可视化）
sp_psi_mean <- array(NA_real_,
                      dim = c(n_draws, n_sp_use, n_periods),
                      dimnames = list(draw = NULL, species = species, period = NULL))

message(sprintf("[10a] Computing species trends: %d draws × %d species", n_draws, n_sp_use))

for (d in seq_len(n_draws)) {
  for (sp in seq_len(n_sp_use)) {
    # 该物种在所有 site 的 psi 时间序列
    if (length(psi_dim) >= 4) {
      psi_sp <- psi_samples[draw_idx[d], sp, , ]  # sites × periods
    } else {
      psi_sp <- psi_samples[draw_idx[d], sp, ]     # sites（单 period）
      psi_sp <- matrix(psi_sp, ncol = 1)
    }

    # 跨 site 取均值（每 period 的空间均值占有率）
    psi_mean_t <- colMeans(psi_sp, na.rm = TRUE)
    sp_psi_mean[d, sp, ] <- psi_mean_t

    if (sum(!is.na(psi_mean_t)) >= 3) {
      # Theil-Sen 斜率
      sp_trend_arr[d, sp, "theil_sen"] <- theil_sen_slope(psi_mean_t)
      # OLS 斜率
      sp_trend_arr[d, sp, "ols"] <- coef(lm(psi_mean_t ~ seq_len(n_periods)))[2]
    }
  }
  if (d %% 50 == 0) message(sprintf("  [10a] ... draw %d/%d", d, n_draws))
}

# ── 10b. 汇总物种趋势后验 ──────────────────────────────────────────
sp_trend_summary <- list()
for (method in c("theil_sen", "ols")) {
  for (sp in seq_len(n_sp_use)) {
    vals <- sp_trend_arr[, sp, method]
    vals <- vals[!is.na(vals)]
    if (length(vals) < 10) next

    # 后验统计
    q <- quantile(vals, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)

    # 趋势方向后验概率：P(slope > 0) 和 P(slope < 0)
    p_pos <- mean(vals > 0)
    p_neg <- mean(vals < 0)

    sp_trend_summary[[length(sp_trend_summary) + 1]] <- tibble(
      species = species[sp],
      method  = method,
      mean    = mean(vals),
      sd      = sd(vals),
      q025    = q[1], q25 = q[2], median = q[3], q75 = q[4], q975 = q[5],
      p_positive = mean(vals > 0),
      p_negative = mean(vals < 0)
    )
  }
}
sp_trend_df <- bind_rows(sp_trend_summary)
write_csv(sp_trend_df, v3_file("results", paste0("table_species_trend_", run_label)))

# ── 10c. 物种分组：扩张 / 稳定 / 收缩 ──────────────────────────────
# 使用 Theil-Sen 结果（更稳健），基于后验概率分类
# 扩张：P(slope > 0) > 0.95
# 收缩：P(slope < 0) > 0.95
# 稳定：其余

sp_classify <- sp_trend_df |>
  filter(method == "theil_sen") |>
  mutate(
    trend_class = case_when(
      p_positive > 0.95 ~ "expanding",
      p_negative > 0.95 ~ "contracting",
      TRUE              ~ "stable"
    ),
    # 更细致的分级
    trend_grade = case_when(
      p_positive > 0.99 ~ "strong_expanding",
      p_positive > 0.95 ~ "expanding",
      p_negative > 0.99 ~ "strong_contracting",
      p_negative > 0.95 ~ "contracting",
      TRUE              ~ "stable"
    )
  )

# Mann-Kendall 检验（基于后验均值时间序列，补充后验概率）
mk_species <- list()
for (sp in seq_len(n_sp_use)) {
  psi_mean_ts <- apply(sp_psi_mean[, sp, ], 2, mean, na.rm = TRUE)
  if (sum(!is.na(psi_mean_ts)) >= 4) {
    mk <- mann_kendall_test(psi_mean_ts)
    mk_species[[length(mk_species) + 1]] <- tibble(
      species = species[sp],
      mk_S = mk$S, mk_tau = mk$tau, mk_p = mk$p_value
    )
  }
}
mk_species_df <- bind_rows(mk_species)

# 合并分类 + MK 检验
sp_classify <- sp_classify |>
  left_join(mk_species_df %>% mutate(species = as.character(species)), by = "species")

write_csv(sp_classify, v3_file("results", paste0("table_species_trend_classify_", run_label)))

# 分组统计
class_table <- sp_classify |>
  count(trend_class, name = "n_species") |>
  mutate(pct = round(n_species / sum(n_species) * 100, 1))
message(sprintf("[10c] Species classification: %s",
  paste(class_table$trend_class, "=", class_table$n_species, collapse = "; ")))

# ── 10d. 物种趋势 × 性状关联（Q4 核心）──────────────────────────────
# 扩张/收缩物种有什么性状特征？

sp_trait_trend <- sp_classify |>
  select(species, trend_class, trend_grade, mean, p_positive, p_negative, mk_tau, mk_p) |>
  rename(trend_slope = mean) |>
  left_join(trait_ext %>% mutate(species = as.character(species)), by = "species")

write_csv(sp_trait_trend, v3_file("results", paste0("table_species_trend_traits_", run_label)))

# 性状与趋势的相关性（Spearman，后验传播）
trait_corr_results <- list()
for (d in seq_len(min(n_draws, 100))) {  # 100 draws 传播不确定性
  for (trait_name in all_trait_vars) {
    if (!trait_name %in% names(trait_ext)) next
    trait_vals <- trait_ext[[trait_name]]
    if (all(is.na(trait_vals))) next

    # 当次 draw 的物种趋势
    slopes <- sp_trend_arr[d, , "theil_sen"]
    names(slopes) <- species

    # 对齐
    common <- intersect(names(slopes), trait_ext$species[!is.na(trait_vals)])
    if (length(common) < 10) next

    slopes_c <- slopes[common]
    traits_c <- trait_ext[[trait_name]][match(common, trait_ext$species)]

    rho <- cor(slopes_c, traits_c, use = "complete.obs", method = "spearman")

    trait_corr_results[[length(trait_corr_results) + 1]] <- tibble(
      draw = d, trait = trait_name, spearman_rho = rho
    )
  }
}
trait_corr_df <- bind_rows(trait_corr_results)
trait_corr_summary <- trait_corr_df |>
  group_by(trait) |>
  summarise(
    rho_mean = mean(spearman_rho, na.rm = TRUE),
    rho_q025 = quantile(spearman_rho, 0.025, na.rm = TRUE),
    rho_q975 = quantile(spearman_rho, 0.975, na.rm = TRUE),
    p_positive = mean(spearman_rho > 0, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(trait_corr_summary, v3_file("results", paste0("table_trait_trend_correlation_", run_label)))
message(sprintf("[10d] Trait-trend correlations: %s",
  paste(trait_corr_summary$trait, "=",
        round(trait_corr_summary$rho_mean, 3), collapse = ", ")))

# ── 10e. 物种趋势 × 环境驱动因子 ────────────────────────────────────
# 扩张物种集中在哪里？什么环境条件？

if (!is.null(grid_env) && nrow(grid_env) > 0) {
  # 计算每物种在扩张/收缩网格中的占有率
  # 策略：用最终 period (P5) 的占有率 × 网格环境变量
  # 物种趋势 vs 其分布区的环境均值

  sp_env_trend <- list()
  # 使用后验均值的最终 period
  psi_final_mean <- apply(psi_samples[draw_idx, , , n_periods], c(2, 3), mean, na.rm = TRUE)
  # psi_final_mean: species × sites

  # 每物种的环境加权均值（按 psi 加权）
  env_driver_vars <- intersect(
    c("bio4", "bio7", "bio11", "bio13", "elev_mean", "hfi_mean",
      "landcover_built", "landcover_cropland", "centroid_lat", "centroid_lon"),
    names(grid_env)
  )

  for (sp in seq_len(min(n_sp_use, 200))) {  # 限制物种数
    psi_sp <- psi_final_mean[sp, ]
    if (all(is.na(psi_sp)) || sum(psi_sp > 0.01) < 5) next

    env_weighted <- numeric(length(env_driver_vars))
    names(env_weighted) <- env_driver_vars
    for (v in env_driver_vars) {
      if (v %in% names(grid_env)) {
        env_vals <- grid_env[[v]]
        env_weighted[v] <- sum(psi_sp * env_vals, na.rm = TRUE) / sum(psi_sp, na.rm = TRUE)
      }
    }

    sp_env_trend[[length(sp_env_trend) + 1]] <- tibble(
      species = species[sp],
      !!!env_weighted
    )
  }
  sp_env_df <- bind_rows(sp_env_trend)

  # 合并物种趋势
  sp_env_trend_df <- sp_classify |>
    select(species, trend_class, trend_slope = mean) |>
    left_join(sp_env_df %>% mutate(species = as.character(species)), by = "species")

  write_csv(sp_env_trend_df, v3_file("results", paste0("table_species_env_trend_", run_label)))

  # 物种趋势 vs 环境的 Spearman 相关
  env_corr_results <- list()
  for (v in env_driver_vars) {
    if (!v %in% names(sp_env_trend_df)) next
    vals <- sp_env_trend_df[[v]]
    slopes <- sp_env_trend_df$trend_slope
    ok <- !is.na(vals) & !is.na(slopes)
    if (sum(ok) < 10) next

    rho <- cor(slopes[ok], vals[ok], method = "spearman")
    env_corr_results[[length(env_corr_results) + 1]] <- tibble(
      env_var = v, spearman_rho = rho, n_species = sum(ok)
    )
  }
  env_corr_df <- bind_rows(env_corr_results)
  write_csv(env_corr_df, v3_file("results", paste0("table_env_trend_correlation_", run_label)))
  message(sprintf("[10e] Environment-trend correlations: %s",
    paste(env_corr_df$env_var, "=", round(env_corr_df$spearman_rho, 3), collapse = ", ")))
} else {
  message("[10e] Skipping environment-trend analysis (grid_env not available)")
}

# ── 10f. 扩张/收缩物种的热点地图数据 ────────────────────────────────
# 对每个网格：扩张物种数 vs 收缩物种数 → 净趋势指数

expanding_sp  <- sp_classify |> filter(trend_class == "expanding")  |> pull(species)
contracting_sp <- sp_classify |> filter(trend_class == "contracting") |> pull(species)

if (length(expanding_sp) > 0 || length(contracting_sp) > 0) {
  # 使用后验均值的 psi
  psi_mean_arr <- apply(psi_samples[draw_idx, , , ], c(2, 3, 4), mean, na.rm = TRUE)
  # psi_mean_arr: species × sites × periods

  # 每个网格每 period 的扩张/收缩物种占有率之和
  hotspot_df <- list()
  for (t in seq_len(n_periods)) {
    # 扩张物种占有率之和
    if (length(expanding_sp) > 0) {
      exp_idx <- which(species %in% expanding_sp)
      exp_richness <- colSums(psi_mean_arr[exp_idx, , t], na.rm = TRUE)
    } else {
      exp_richness <- rep(0, n_sites)
    }

    # 收缩物种占有率之和
    if (length(contracting_sp) > 0) {
      con_idx <- which(species %in% contracting_sp)
      con_richness <- colSums(psi_mean_arr[con_idx, , t], na.rm = TRUE)
    } else {
      con_richness <- rep(0, n_sites)
    }

    hotspot_df[[length(hotspot_df) + 1]] <- tibble(
      grid_cell = sites,
      period = paste0("P", t),
      expanding_richness = exp_richness,
      contracting_richness = con_richness,
      net_trend_index = exp_richness - con_richness,
      turnover_balance = ifelse(exp_richness + contracting_richness > 0,
                                 (exp_richness - contracting_richness) / (exp_richness + contracting_richness),
                                 0)
    )
  }
  hotspot_all <- bind_rows(hotspot_df)
  write_csv(hotspot_all, v3_file("results", paste0("table_species_hotspot_", run_label)))
  message(sprintf("[10f] Hotspot map data: %d grids × %d periods", n_sites, n_periods))
}

message(sprintf("[10] Species-level analysis complete: %d expanding, %d stable, %d contracting",
  sum(sp_classify$trend_class == "expanding"),
  sum(sp_classify$trend_class == "stable"),
  sum(sp_classify$trend_class == "contracting")))

log_time("05", "DONE")
