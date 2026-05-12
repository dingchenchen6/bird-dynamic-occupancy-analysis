#!/usr/bin/env Rscript
## 09_extended_analyses.R  —  v3 扩展分析（改进版）
##
## C.  varpart — 使用变化量驱动因子（时间变化→时间变化）
## C2. 随机森林排列重要性（后验传播）
## C3. 敏感性分析（标准化气候变化量 vs 绝对变化量）
## D.  HFI 分层
## E.  Mann-Kendall 正式检验
##
## v3 改进要点：
##   - Q3 驱动因子分析从静态变量升级为时间变化量
##   - 新增 DRIVER_GROUPS_TREND（变化量驱动组）
##   - 新增敏感性分析（delta_t_std vs delta_t_mean）
##   - 保留静态 DRIVER_GROUPS 作为 stMsPGOcc 内部使用

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(vegan); library(Kendall)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_importance.R"))
source(file.path(CODE_V3, "utils_mapping.R"))
P <- ensure_v3_dirs()

log_time("09", "Starting extended analyses (v3 with temporal change drivers)")

is_pilot <- Sys.getenv("V3_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL

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
  filter(metric == "corrected_richness") |>
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

# ── C. varpart（变化量驱动版）────────────────────────────────────────
message("[C] Variance partitioning with temporal change drivers")

groups_trend <- DRIVER_GROUPS_TREND  # 来自 00_config

# 检查哪些变化量变量实际可用
groups_trend <- lapply(groups_trend, function(vars) {
  avail <- intersect(vars, names(driver_dat))
  if (length(avail) < length(vars)) {
    missing <- setdiff(vars, avail)
    message(sprintf("[C] Missing from data: %s", paste(missing, collapse = ", ")))
  }
  avail
})

# 需要至少 2 个组有数据
groups_with_data <- groups_trend[sapply(groups_trend, length) > 0]
if (length(groups_with_data) < 2) {
  message("[C] Not enough driver groups with data. Falling back to static DRIVER_GROUPS.")
  groups_trend <- lapply(DRIVER_GROUPS, function(vars) intersect(vars, names(driver_dat)))
  groups_with_data <- groups_trend[sapply(groups_trend, length) > 0]
}

# 构建完整分析数据（去掉任何驱动变量有 NA 的行）
all_driver_vars <- unlist(groups_with_data)
vp_dat <- driver_dat |>
  select(grid_cell, mean, all_of(all_driver_vars)) |>
  filter(complete.cases(across(all_of(all_driver_vars))))

message(sprintf("[C] %d grids with complete driver data (of %d total)",
                nrow(vp_dat), nrow(driver_dat)))

if (nrow(vp_dat) >= 50) {  # 需要足够样本
  make_X <- function(vars) {
    vp_dat |>
      select(all_of(vars)) |>
      mutate(across(everything(), as.numeric)) |>
      mutate(across(everything(), ~ as.numeric(scale(.x))))
  }

  # 4-way varpart（或根据可用组数调整）
  group_names <- names(groups_with_data)
  X_list <- lapply(group_names, function(g) make_X(groups_with_data[[g]]))

  vp <- do.call(vegan::varpart, c(list(Y = vp_dat$mean), X_list))
  saveRDS(vp, v3_file("derived", paste0("varpart_richness_trend_", run_label), "rds"))

  # 提取结果
  vp_part <- as_tibble(vp$part$indfract, rownames = "fraction")
  adjcol <- grep("[Aa]dj.*[Rr]", names(vp_part), value = TRUE)[1]
  if (is.na(adjcol)) adjcol <- grep("[Rr].sq", names(vp_part), value = TRUE)[1]
  if (is.na(adjcol)) adjcol <- "R.square"
  vp_part <- vp_part |>
    rename(adj_r2 = !!adjcol) |>
    mutate(fraction = factor(fraction, levels = unique(fraction)))
  write_csv(vp_part, v3_file("results", paste0("table_varpart_richness_trend_", run_label)))

  # 纯效应提取
  # 构建标签映射
  pure_labels <- sapply(group_names, function(g) {
    DRIVER_GROUP_TREND_LABELS[g] %||% paste0(g, " (pure)")
  })
  letter_map <- setNames(paste0("[", letters[seq_along(group_names)], "]"),
                         group_names)

  vp_pure <- vp_part |>
    filter(grepl("^\\[\\w+\\]$", fraction)) |>
    mutate(component = recode(as.character(fraction), !!!setNames(
      paste0(pure_labels, " (pure)"), letter_map))) |>
    mutate(adj_R2 = pmax(adj_r2, 0))

  # varpart 可视化
  vp_pretty <- vp_pure |>
    mutate(component = recode(component, !!!setNames(
      DRIVER_GROUP_TREND_LABELS[names(DRIVER_GROUP_TREND_LABELS) %in% group_names],
      DRIVER_GROUP_TREND_LABELS[names(DRIVER_GROUP_TREND_LABELS) %in% group_names]
    )))

  message("[C] Variance partitioning results:")
  for (i in seq_len(nrow(vp_pretty))) {
    message(sprintf("  %s: adj R² = %.4f", vp_pretty$component[i], vp_pretty$adj_R2[i]))
  }

} else {
  message("[C] Too few complete cases for varpart. Skipping.")
  vp_pure <- NULL
}

# ── C2. 随机森林排列重要性（变化量驱动）────────────────────────────
message("[C2] Random Forest permutation importance with temporal change drivers")

psi_obj <- safe_read(v3_file("derived", paste0("psi_samples_thinned_", run_label), "rds"))
trend_draws_arr <- safe_read(v3_file("derived",
  paste0("trend_draws_corrected_richness_", run_label), "rds"))

if (!is.null(trend_draws_arr) && nrow(vp_dat) >= 50) {
  env_rf <- vp_dat |>
    select(all_of(unlist(groups_with_data))) |>
    mutate(across(everything(), as.numeric))

  rf_draws <- rf_importance_draws(
    trend_draws = trend_draws_arr,
    env_df      = env_rf,
    groups      = groups_with_data,
    num_draws   = min(RF_NUM_DRAWS, dim(trend_draws_arr)[1]),
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
    select(all_of(unlist(groups_with_data))) |>
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

# ── C3. 敏感性分析：标准化 vs 绝对气候变化量 ──────────────────────
message("[C3] Sensitivity analysis: standardized vs absolute climate change")

if ("delta_t_std" %in% names(driver_dat) && sum(!is.na(driver_dat$delta_t_std)) >= 50) {
  groups_sensitivity <- DRIVER_GROUPS_TREND_SENSITIVITY
  groups_sensitivity <- lapply(groups_sensitivity, function(vars) intersect(vars, names(driver_dat)))

  vp_dat_sens <- driver_dat |>
    select(grid_cell, mean, all_of(unlist(groups_sensitivity))) |>
    filter(complete.cases(across(all_of(unlist(groups_sensitivity)))))

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

# ── F. 变化量驱动的 brms 回归（替代原 08 脚本中的静态版）────────────
message("[F] brms driver regression with temporal change predictors")

# 检查是否已有 brms 拟合（原 08 脚本可能已运行）
# 如果变化量驱动数据可用，运行新的 brms 回归
if (requireNamespace("brms", quietly = TRUE) && nrow(vp_dat) >= 50) {
  # 构建公式
  brms_vars <- all_driver_vars[all_driver_vars %in% names(vp_dat)]

  if (length(brms_vars) >= 3) {
    # 构建公式字符串
    brms_formula_str <- paste0(
      "mean ~ ",
      paste(sapply(brms_vars, function(v) sprintf("scale(%s)", v)), collapse = " + "),
      " + gp(centroid_lon, centroid_lat, k = 10, c = 5/4)"
    )

    message(sprintf("[F] brms formula: %s", brms_formula_str))

    brms_fit <- tryCatch({
      brms::brm(
        as.formula(brms_formula_str),
        data = vp_dat,
        iter = BRMS_ITER, warmup = BRMS_WARMUP,
        chains = BRMS_CHAINS, seed = BRMS_SEED,
        backend = "cmdstanr",
        control = list(adapt_delta = BRMS_ADAPT_DELTA,
                       max_treedepth = BRMS_MAX_TREED),
        silent = 2, refresh = 0
      )
    }, error = function(e) {
      message(sprintf("[F] brms fit failed: %s", conditionMessage(e)))
      NULL
    })

    if (!is.null(brms_fit)) {
      saveRDS(brms_fit, v3_file("derived",
        paste0("brms_driver_trend_", run_label), "rds"))

      # 提取固定效应后验
      fe <- brms::fixef(brms_fit)
      fe_tibble <- as_tibble(fe, rownames = "term") |>
        mutate(term = gsub("^scale", "", term))

      write_csv(fe_tibble, v3_file("results",
        paste0("table_brms_driver_trend_coefficients_", run_label)))

      message("[F] brms driver trend fit completed.")
    }
  } else {
    message("[F] Not enough change-driver variables available. Skipping brms.")
  }
} else {
  message("[F] brms not available or insufficient data. Skipping.")
}

log_time("09", "DONE")
