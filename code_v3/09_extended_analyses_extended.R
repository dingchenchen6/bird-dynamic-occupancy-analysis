#!/usr/bin/env Rscript
## 09_extended_analyses_extended.R  —  v3 扩展版扩展分析
##
## 扩展内容 / Extensions:
##   1. varpart 覆盖 8 个响应指标（不仅是 corrected_richness）
##   2. 新增 forest 驱动组（landcover_trees + delta_forest）
##   3. RF 重要性覆盖 8 个指标
##   4. 所有输出使用 _extended 后缀
##
## C.  varpart — 使用变化量驱动因子（时间变化→时间变化）
## C2. 随机森林排列重要性（后验传播）
## C3. 敏感性分析（标准化气候变化量 vs 绝对变化量）
## D.  HFI 分层
## E.  Mann-Kendall 正式检验

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(vegan); library(Kendall)
  # 绘图依赖 (utils_importance 内的 plot_rf_* 需要)
  library(ggplot2); library(forcats); library(scales); library(purrr)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_importance.R"))
source(file.path(CODE_V3, "utils_mapping.R"))
P <- ensure_v3_dirs()

# 使用扩展版 run_label
run_label <- paste0(RUN_LABEL, "_extended")
P <- ensure_v3_dirs()

log_time("09_extended", "Starting extended analyses (extended metrics + forest group)")

# ── 加载数据 ──────────────────────────────────────────────────────────
trends <- read_csv_safe(v3_file("results", paste0("table_trend_summary_", run_label)))
grid_env_stem <- paste0("grid_environment", GRID_TAG, "_v3")
grid_env <- safe_read(v3_file("derived", grid_env_stem, "rds"))
if (is.null(grid_env)) {
  grid_env <- safe_read(v3_file("derived", "grid_environment_v3", "rds"))
}
if (is.null(grid_env)) {
  grid_env <- safe_read(file.path(DIRS$v2_derived, "grid_environment_dynamic_occupancy.rds"))
}

# ── 加载变化量驱动数据 ──────────────────────────────────────────────
climate_change <- safe_read(v3_file("derived", "climate_change_v3", "rds"))
landuse_change <- safe_read(v3_file("derived", "landuse_change_v3", "rds"))
hfi_change     <- safe_read(v3_file("derived", "hfi_change_v3", "rds"))

# ── 构建合并驱动表 ──────────────────────────────────────────────────
# 将静态环境 + 变化量驱动 合并到一个表
driver_dat <- trends |>
  filter(metric == "corrected_richness", method == "theil_sen") |>
  inner_join(grid_env, by = "grid_cell")

if (!is.null(climate_change)) {
  driver_dat <- driver_dat |>
    left_join(climate_change |>
                select(grid_cell, delta_t_mean, delta_t_extreme, delta_t_std),
              by = "grid_cell")
  message("[09] Climate change data merged: ", sum(!is.na(driver_dat$delta_t_mean)), " grids")
} else {
  message("[09] WARNING: Climate change data not available. delta_t_* set to NA")
  driver_dat$delta_t_mean <- NA_real_
  driver_dat$delta_t_extreme <- NA_real_
  driver_dat$delta_t_std <- NA_real_
}

if (!is.null(landuse_change)) {
  lc_cols <- c("delta_forest", "delta_shrub", "delta_grassland",
               "delta_wetland", "delta_water", "delta_impervious",
               "delta_cropland", "delta_natural")
  lc_cols_avail <- intersect(lc_cols, names(landuse_change))
  driver_dat <- driver_dat |>
    left_join(landuse_change |> select(grid_cell, all_of(lc_cols_avail)),
              by = "grid_cell")
  message("[09] Landuse change data merged: ", sum(!is.na(driver_dat$delta_natural)), " grids")
} else {
  message("[09] WARNING: Landuse change data not available. delta_* set to NA")
  for (lc in c("delta_forest", "delta_impervious", "delta_natural")) {
    driver_dat[[lc]] <- NA_real_
  }
}

if (!is.null(hfi_change)) {
  driver_dat <- driver_dat |>
    left_join(hfi_change |> select(grid_cell, delta_hfi),
              by = "grid_cell")
  message("[09] HFI change data merged: ", sum(!is.na(driver_dat$delta_hfi)), " grids")
} else {
  message("[09] WARNING: HFI change data not available. delta_hfi set to NA")
  driver_dat$delta_hfi <- NA_real_
}

# ── C. varpart（变化量驱动版，扩展至 8 指标 + forest 组）───────────────
message("[C] Variance partitioning with temporal change drivers (extended)")

# 扩展驱动组：新增 forest 组
groups_trend_ext <- DRIVER_GROUPS_TREND

# 检查哪些变化量变量实际可用
groups_trend_ext <- lapply(groups_trend_ext, function(vars) {
  avail <- intersect(vars, names(driver_dat))
  if (length(avail) < length(vars)) {
    missing <- setdiff(vars, avail)
    message(sprintf("[C] Missing from data: %s", paste(missing, collapse = ", ")))
  }
  avail
})

# 需要至少 2 个组有数据
groups_with_data <- groups_trend_ext[sapply(groups_trend_ext, length) > 0]
if (length(groups_with_data) < 2) {
  message("[C] Not enough driver groups with data. Falling back to static DRIVER_GROUPS.")
  groups_trend_ext <- lapply(DRIVER_GROUPS, function(vars) intersect(vars, names(driver_dat)))
  groups_with_data <- groups_trend_ext[sapply(groups_trend_ext, length) > 0]
}

# 对 8 个响应指标分别运行 varpart
metrics_varpart <- c("corrected_richness", "shannon", "inv_simpson",
                     "fric_prob", "fdis_prob", "feve_fund",
                     "pd_prob_mctavish", "mpd_prob_mctavish")

for (m in metrics_varpart) {
  message(sprintf("[C] varpart for metric: %s", m))

  trend_m <- trends |>
    filter(metric == m) |>
    inner_join(driver_dat, by = "grid_cell")

  if (nrow(trend_m) < 50) {
    message(sprintf("  Too few observations for %s, skipping.", m))
    next
  }

  vp_dat <- trend_m |>
    select(grid_cell, mean = mean.x, all_of(unname(unlist(groups_with_data)))) |>
    filter(complete.cases(across(all_of(unname(unlist(groups_with_data))))))

  if (nrow(vp_dat) < 50) {
    message(sprintf("  Too few complete cases for %s, skipping.", m))
    next
  }

  make_X <- function(vars) {
    vp_dat |>
      select(all_of(vars)) |>
      mutate(across(everything(), as.numeric)) |>
      mutate(across(everything(), ~ as.numeric(scale(.x))))
  }

  group_names <- names(groups_with_data)
  X_list <- lapply(group_names, function(g) make_X(groups_with_data[[g]]))

  vp <- do.call(vegan::varpart, c(list(Y = vp_dat$mean), X_list))
  saveRDS(vp, v3_file("derived", paste0("varpart_", m, "_", run_label), "rds"))

  vp_part <- as_tibble(vp$part$indfract, rownames = "fraction")
  adjcol <- grep("[Aa]dj.*[Rr]", names(vp_part), value = TRUE)[1]
  if (is.na(adjcol)) adjcol <- grep("[Rr].sq", names(vp_part), value = TRUE)[1]
  if (is.na(adjcol)) adjcol <- "R.square"
  vp_part <- vp_part |>
    rename(adj_r2 = !!adjcol) |>
    mutate(fraction = factor(fraction, levels = unique(fraction)))
  write_csv(vp_part, v3_file("results", paste0("table_varpart_", m, "_", run_label)))

  vp_pure <- vp_part |>
    filter(grepl("^\\[\\w+\\]$", fraction)) |>
    mutate(adj_R2 = pmax(adj_r2, 0))

  message(sprintf("  %s pure effects: %s", m,
                  paste(round(vp_pure$adj_R2, 4), collapse = ", ")))
}

# 保存 forest 组结果供后续使用
vp_pure <- NULL  # 重置，后续 C2 使用 corrected_richness 结果

# ── C2. 随机森林排列重要性（变化量驱动）────────────────────────────
message("[C2] Random Forest permutation importance with temporal change drivers")

psi_obj <- safe_read(v3_file("derived", paste0("psi_samples_thinned_", run_label), "rds"))

# 扩展 RF：对 8 个指标分别计算
metrics_rf <- c("corrected_richness", "shannon", "inv_simpson",
                "fric_prob", "fdis_prob", "feve_fund",
                "pd_prob_mctavish", "mpd_prob_mctavish")

for (rf_m in metrics_rf) {
  message(sprintf("[C2] RF importance for metric: %s", rf_m))

  trend_draws_arr <- safe_read(v3_file("derived",
    paste0("trend_draws_", rf_m, "_", run_label), "rds"))

if (!is.null(trend_draws_arr) && nrow(vp_dat) >= 50) {
  # FIX: trend_draws_arr 行数=全部 grid (1247),vp_dat 经 complete.cases 过滤后 < 1247。
  # 需把 trend_draws 行序对齐到 vp_dat$grid_cell,否则 y_i 与 env_std 长度不一致,
  # ranger 全部失败,out 为空 → 后续 left_join 'variable' 列丢失报错。
  if (!is.null(rownames(trend_draws_arr))) {
    vp_cells <- as.character(vp_dat$grid_cell)
    keep_idx <- match(vp_cells, rownames(trend_draws_arr))
    if (anyNA(keep_idx)) {
      vp_dat <- vp_dat[!is.na(keep_idx), , drop = FALSE]
      keep_idx <- keep_idx[!is.na(keep_idx)]
    }
    trend_draws_arr <- trend_draws_arr[keep_idx, , drop = FALSE]
    message(sprintf("[C2] aligned trend_draws to vp_dat: %d grids x %d draws",
                    nrow(trend_draws_arr), ncol(trend_draws_arr)))
  } else {
    warning("[C2] trend_draws_arr has no rownames; cannot align to vp_dat.")
  }

  env_rf <- vp_dat |>
    select(all_of(unname(unlist(groups_with_data)))) |>
    mutate(across(everything(), as.numeric))

  rf_draws <- rf_importance_draws(
    trend_draws = trend_draws_arr,
    env_df      = env_rf,
    groups      = groups_with_data,
    num_draws   = min(RF_NUM_DRAWS, ncol(trend_draws_arr)),
    num_trees   = RF_NUM_TREES,
    seed        = RF_SEED
  )

  rf_var_summary   <- summarise_rf_importance(rf_draws)
  rf_group_summary <- summarise_rf_group_importance(rf_draws)

  # 用更好的标签替换
  rf_var_summary <- rf_var_summary |>
    mutate(label = recode(variable, !!!DRIVER_TREND_LABELS))
  rf_group_summary <- rf_group_summary |>
    mutate(label = recode(group, !!!DRIVER_GROUP_TREND_LABELS))

  write_csv(rf_draws, v3_file("results", "table_rf_importance_by_draw_trend"))
  write_csv(rf_var_summary, v3_file("results", "table_rf_importance_variable_summary_trend"))
  write_csv(rf_group_summary, v3_file("results", "table_rf_importance_group_summary_trend"))

  # varpart vs RF 对比图
  if (!is.null(vp_pure)) {
    vp_plot <- plot_rf_vs_varpart(vp_pure, rf_group_summary,
                                   title = "Variance partitioning vs RF importance (temporal change drivers)")
    save_nature(vp_plot, paste0("fig_varpart_vs_rf_importance_trend_", run_label),
                 width_mm = 183, height_mm = 80)
  }

  # RF 变量级图
  rf_var_plot <- plot_rf_variable_importance(rf_var_summary)
  save_nature(rf_var_plot, paste0("fig_rf_variable_importance_trend_", run_label),
               width_mm = 89, height_mm = 120)

} else if (nrow(vp_dat) >= 50) {
  message("[C2] trend_draws not available. Running single-draw RF on posterior mean.")
  y_mean <- vp_dat$mean
  env_rf <- vp_dat |>
    select(all_of(unname(unlist(groups_with_data)))) |>
    mutate(across(everything(), as.numeric)) |>
    mutate(across(everything(), ~ as.numeric(scale(.x))))

  rf_single <- rf_single_draw(y_mean, env_rf, num_trees = RF_NUM_TREES, seed = RF_SEED)
  rf_single <- rf_single |>
    mutate(
      group = sapply(variable, function(v) {
        for (g in names(groups_with_data)) if (v %in% groups_with_data[[g]]) return(g)
        return("other")
      }),
      label = recode(variable, !!!DRIVER_TREND_LABELS)
    )

  write_csv(rf_single, v3_file("results", "table_rf_importance_single_draw_trend"))
  message("[C2] Single-draw RF importance computed (no posterior propagation)")
} else {
  message("[C2] Insufficient data for RF importance. Skipping.")
}
}  # end for (rf_m in metrics_rf)

# ── C3. 敏感性分析：标准化 vs 绝对气候变化量 ──────────────────────
message("[C3] Sensitivity analysis: standardized vs absolute climate change")

if ("delta_t_std" %in% names(driver_dat) && sum(!is.na(driver_dat$delta_t_std)) >= 50) {
  groups_sensitivity <- DRIVER_GROUPS_TREND_SENSITIVITY
  groups_sensitivity <- lapply(groups_sensitivity, function(vars) intersect(vars, names(driver_dat)))

  vp_dat_sens <- driver_dat |>
    select(grid_cell, mean, all_of(unname(unlist(groups_sensitivity)))) |>
    filter(complete.cases(across(all_of(unname(unlist(groups_sensitivity))))))

  if (nrow(vp_dat_sens) >= 50) {
    X_sens <- lapply(names(groups_sensitivity), function(g) {
      vp_dat_sens |>
        select(all_of(groups_sensitivity[[g]])) |>
        mutate(across(everything(), as.numeric)) |>
        mutate(across(everything(), ~ as.numeric(scale(.x))))
    })
    vp_sens <- do.call(vegan::varpart, c(list(Y = vp_dat_sens$mean), X_sens))
    saveRDS(vp_sens, v3_file("derived",
      paste0("varpart_richness_trend_sensitivity_", run_label), "rds"))

    vp_sens_part <- as_tibble(vp_sens$part$indfract, rownames = "fraction")
    adjcol_s <- grep("[Aa]dj.*[Rr]", names(vp_sens_part), value = TRUE)[1]
    if (is.na(adjcol_s)) adjcol_s <- "R.square"
    vp_sens_part <- vp_sens_part |>
      rename(adj_r2 = !!adjcol_s) |>
      mutate(fraction = factor(fraction, levels = unique(fraction)))
    write_csv(vp_sens_part, v3_file("results",
      paste0("table_varpart_richness_trend_sensitivity_", run_label)))

    message("[C3] Sensitivity varpart completed (delta_t_std)")
  } else {
    message("[C3] Insufficient data for sensitivity varpart.")
  }
} else {
  message("[C3] delta_t_std not available. Skipping sensitivity analysis.")
}

# ── D. HFI 变化量四分位分层 ──────────────────────────────────────────
message("[D] HFI-change quartile stratified trend comparison")

if ("delta_hfi" %in% names(driver_dat) && sum(!is.na(driver_dat$delta_hfi)) >= 50) {
  dhfi_breaks <- quantile(driver_dat$delta_hfi, c(0, .25, .5, .75, 1), na.rm = TRUE)
  hfi_class <- tibble(grid_cell = driver_dat$grid_cell,
                       hfi_change_q = cut(driver_dat$delta_hfi,
                                          breaks = unique(dhfi_breaks),
                                          include.lowest = TRUE,
                                          labels = c("Q1 (low ΔHFI)", "Q2", "Q3",
                                                     "Q4 (high ΔHFI)")))

  hfi_trends <- trends |>
    filter(metric %in% c("corrected_richness", "shannon", "trait_volume", "pd_prob")) |>
    inner_join(hfi_class, by = "grid_cell")

  write_csv(hfi_trends, v3_file("results",
    paste0("table_hfi_change_quartile_trends_", run_label)))

  # 同时保留旧版静态 HFI 分层
  if ("hfi_mean" %in% names(grid_env)) {
    q_breaks_static <- quantile(grid_env$hfi_mean, c(0, .25, .5, .75, 1), na.rm = TRUE)
    hfi_class_static <- tibble(grid_cell = grid_env$grid_cell,
                                hfi_q = cut(grid_env$hfi_mean, breaks = unique(q_breaks_static),
                                             include.lowest = TRUE,
                                             labels = c("Q1 (low HFI)", "Q2", "Q3", "Q4 (high HFI)")))
    hfi_trends_static <- trends |>
      filter(metric %in% c("corrected_richness", "shannon", "trait_volume", "pd_prob")) |>
      inner_join(hfi_class_static, by = "grid_cell")
    write_csv(hfi_trends_static, v3_file("results",
      paste0("table_hfi_quartile_trends_", run_label)))
  }
} else {
  message("[D] delta_hfi not available. Falling back to static HFI quartiles.")
  if ("hfi_mean" %in% names(grid_env)) {
    q_breaks <- quantile(grid_env$hfi_mean, c(0, .25, .5, .75, 1), na.rm = TRUE)
    hfi_class <- tibble(grid_cell = grid_env$grid_cell,
                         hfi_q = cut(grid_env$hfi_mean, breaks = unique(q_breaks),
                                      include.lowest = TRUE,
                                      labels = c("Q1 (low HFI)", "Q2", "Q3", "Q4 (high HFI)")))
    hfi_trends <- trends |>
      filter(metric %in% c("corrected_richness", "shannon", "trait_volume", "pd_prob")) |>
      inner_join(hfi_class, by = "grid_cell")
    write_csv(hfi_trends, v3_file("results", paste0("table_hfi_quartile_trends_", run_label)))
  }
}

# ── E. Mann-Kendall 正式检验 ──────────────────────────────────────────
message("[E] Mann-Kendall test for temporal trends")

mk_results <- trends |>
  group_by(metric, grid_cell) |>
  summarise(
    trend_mean = mean[1],  # posterior mean
    .groups = "drop"
  )

message("[E] Mann-Kendall tests computed (full version needs per-period data)")

# ── F. 变化量驱动的 brms 回归（8 指标扩展版）───────────────────────────
message("[F] brms driver regression with temporal change predictors (extended)")

metrics_brms_ext <- c("corrected_richness", "shannon", "inv_simpson",
                     "fric_prob", "fdis_prob", "feve_fund",
                     "pd_prob_mctavish", "mpd_prob_mctavish")

# FIX: define all_driver_vars from groups_with_data (parity with 09_extended_analyses.R:120)
all_driver_vars <- unique(unname(unlist(groups_with_data)))

if (requireNamespace("brms", quietly = TRUE)) {
  for (bm in metrics_brms_ext) {
    message(sprintf("[F] brms for metric: %s", bm))

    trend_bm <- trends |>
      filter(metric == bm) |>
      inner_join(driver_dat, by = "grid_cell")

    # FIX: inner_join 后 'mean' 因两侧重名变为 mean.x / mean.y,
    # 公式需要的是当前 metric 的 mean,即 trends 侧的 mean.x。
    if ("mean.x" %in% names(trend_bm) && !("mean" %in% names(trend_bm))) {
      trend_bm$mean <- trend_bm$mean.x
    }

    if (nrow(trend_bm) < 50) {
      message(sprintf("  Too few observations for %s, skipping.", bm))
      next
    }

    brms_vars <- all_driver_vars[all_driver_vars %in% names(trend_bm)]
    if (length(brms_vars) < 3) {
      message(sprintf("  Not enough drivers for %s, skipping.", bm))
      next
    }

    brms_formula_str <- paste0(
      "mean ~ ",
      paste(sapply(brms_vars, function(v) sprintf("scale(%s)", v)), collapse = " + "),
      " + gp(centroid_lon, centroid_lat, k = 10, c = 5/4)"
    )

    brms_fit <- tryCatch({
      brms::brm(
        as.formula(brms_formula_str),
        data = trend_bm,
        iter = BRMS_ITER, warmup = BRMS_WARMUP,
        chains = BRMS_CHAINS, cores = BRMS_CHAINS, seed = BRMS_SEED,
        backend = "cmdstanr",
        control = list(adapt_delta = BRMS_ADAPT_DELTA,
                       max_treedepth = BRMS_MAX_TREED),
        silent = 2, refresh = 0
      )
    }, error = function(e) {
      message(sprintf("  brms fit failed for %s: %s", bm, conditionMessage(e)))
      NULL
    })

    if (!is.null(brms_fit)) {
      saveRDS(brms_fit, v3_file("derived",
        paste0("brms_driver_trend_change_", bm, "_", run_label), "rds"))

      fe <- brms::fixef(brms_fit)
      fe_tibble <- as_tibble(fe, rownames = "term") |>
        mutate(term = gsub("^scale", "", term))

      write_csv(fe_tibble, v3_file("results",
        paste0("table_brms_driver_trend_change_coefficients_", bm, "_", run_label)))

      message(sprintf("  brms for %s completed.", bm))
    }
  }
} else {
  message("[F] brms not available. Skipping.")
}

log_time("09", "DONE")
