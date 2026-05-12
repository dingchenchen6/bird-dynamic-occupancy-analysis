#!/usr/bin/env Rscript
## utils_importance.R  —  v3 随机森林重要性评估工具
##
## 使用 ranger 进行排列重要性评估，传播后验不确定性
## 与 vegan::varpart 互补：RF 捕捉非线性关系
##
## 参考文献：Morelli et al. 2021 的 diet specialization 方法论启发
## 随机森林方法论：Breiman 2001; ranger: Wright & Ziegler 2017

# ── 加载依赖 ──────────────────────────────────────────────────────────
{
  .root <- Sys.getenv("BIRD_PROJECT_ROOT",
    if (dir.exists(file.path("~", "Documents", "New project",
                             "bird_dynamic_occupancy_analysis")))
      file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis")
    else getwd())
  source(file.path(.root, "code_v3", "utils_paths.R"))
  rm(.root)
}

# ── 随机森林排列重要性（单次后验抽取）─────────────────────────────────
#' 对单个后验抽取拟合 ranger 并返回排列重要性
#'
#' @param y 响应向量（某 draw 的 trend 值）
#' @param X 环境变量 data.frame
#' @param num_trees 树的数量
#' @param seed 随机种子
#' @return tibble: variable, importance
rf_single_draw <- function(y, X, num_trees = 1000, seed = NULL) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("ranger package required. Install with: install.packages('ranger')")
  }

  df <- cbind(trend = y, X)
  fit <- ranger::ranger(
    formula   = trend ~ .,
    data      = df,
    num.trees = num_trees,
    importance = "permutation",
    seed       = seed %||% sample.int(1e6, 1),
    min.node.size = 5,
    verbose   = FALSE
  )

  imp <- ranger::importance(fit)
  tibble(variable = names(imp), importance = as.numeric(imp))
}

# ── 随机森林排列重要性（多后验抽取，传播不确定性）─────────────────────
#' 对多个后验抽取重复 RF，汇总重要性分布
#'
#' @param trend_draws 矩阵 (n_grid × n_draws)：每个 draw 一列 trend 值
#' @param env_df 环境变量 data.frame (n_grid × p)
#' @param groups 变量分组列表（与 varpart 一致）
#' @param num_draws 使用多少后验抽取
#' @param num_trees 每次拟合的树数
#' @param seed 基础随机种子
#' @return tibble: draw, variable, importance, group
rf_importance_draws <- function(trend_draws, env_df, groups = NULL,
                                 num_draws = 100, num_trees = 1000,
                                 seed = 2025) {
  n_draws_avail <- ncol(trend_draws)
  draw_idx <- if (num_draws >= n_draws_avail) {
    seq_len(n_draws_avail)
  } else {
    set.seed(seed)
    sample(n_draws_avail, num_draws)
  }

  # 标准化环境变量
  env_std <- env_df |>
    mutate(across(where(is.numeric), ~ as.numeric(scale(.x))))

  # 分组映射
  var_to_group <- if (!is.null(groups)) {
    purrr::map2_dfr(groups, names(groups), ~ tibble(variable = .x, group = .y))
  } else {
    NULL
  }

  message(sprintf("[rf_importance] Running %d RF fits × %d trees ...",
                  length(draw_idx), num_trees))

  results <- list()
  for (i in seq_along(draw_idx)) {
    d <- draw_idx[i]
    y_i <- trend_draws[, d]
    seed_i <- seed + i

    imp_i <- tryCatch(
      rf_single_draw(y_i, env_std, num_trees = num_trees, seed = seed_i),
      error = function(e) {
        warning(sprintf("RF draw %d failed: %s", d, e$message), call. = FALSE)
        NULL
      }
    )

    if (!is.null(imp_i)) {
      results[[i]] <- imp_i |> mutate(draw = d)
    }

    if (i %% 10 == 0) message(sprintf("  ... %d/%d draws done", i, length(draw_idx)))
  }

  out <- bind_rows(results)

  # 加入分组信息
  if (!is.null(var_to_group)) {
    out <- out |> left_join(var_to_group, by = "variable")
  } else {
    out <- out |> mutate(group = "All")
  }

  out
}

# ── 汇总 RF 重要性 ──────────────────────────────────────────────────
#' 汇总多次 RF 拟合的排列重要性
#'
#' @param rf_draws_df rf_importance_draws() 的输出
#' @param level CRI 水平（默认 0.95）
#' @return tibble: variable, group, mean, median, sd, q_lo, q_hi
summarise_rf_importance <- function(rf_draws_df, level = 0.95) {
  alpha <- 1 - level
  rf_draws_df |>
    group_by(variable, group) |>
    summarise(
      mean   = mean(importance, na.rm = TRUE),
      median = median(importance, na.rm = TRUE),
      sd     = sd(importance, na.rm = TRUE),
      q_lo   = quantile(importance, alpha / 2, na.rm = TRUE),
      q_hi   = quantile(importance, 1 - alpha / 2, na.rm = TRUE),
      n_draws = sum(!is.na(importance)),
      .groups = "drop"
    ) |>
    arrange(desc(mean))
}

# ── 按组聚合重要性 ────────────────────────────────────────────────────
summarise_rf_group_importance <- function(rf_draws_df) {
  # 先按 draw × group 求和
  by_draw_group <- rf_draws_df |>
    group_by(draw, group) |>
    summarise(importance = sum(importance, na.rm = TRUE), .groups = "drop")

  # 再按组汇总
  by_draw_group |>
    group_by(group) |>
    summarise(
      mean   = mean(importance, na.rm = TRUE),
      median = median(importance, na.rm = TRUE),
      sd     = sd(importance, na.rm = TRUE),
      q_lo   = quantile(importance, 0.025, na.rm = TRUE),
      q_hi   = quantile(importance, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(desc(mean))
}

# ── varpart vs RF 对比图（Nature 风格）────────────────────────────────
#' 左面板：varpart 纯分数；右面板：RF 组级重要性（带误差棒）
#'
#' @param varpart_summary tibble: component, adj_R2（来自 09_extended_analyses）
#' @param rf_group_summary summarise_rf_group_importance() 的输出
#' @param title 标题
plot_rf_vs_varpart <- function(varpart_summary, rf_group_summary, title = NULL) {
  # 左面板：varpart
  vp_data <- varpart_summary |>
    mutate(component = recode(component,
      "Climate (pure)" = "Climate",
      "Topo+Habitat (pure)" = "Topo+Habitat",
      "Human (pure)" = "Human",
      "Space (pure)" = "Space"
    )) |>
    rename(mean = adj_R2) |>
    mutate(method = "varpart\n(linear)")

  # 右面板：RF
  rf_data <- rf_group_summary |>
    mutate(method = "RF\n(nonlinear)")

  combined <- bind_rows(vp_data, rf_data)

  ggplot(combined, aes(x = reorder(component, mean), y = mean, fill = method)) +
    geom_col(position = "dodge", width = 0.6, alpha = 0.85) +
    geom_errorbar(aes(ymin = q_lo %||% mean, ymax = q_hi %||% mean),
                  position = position_dodge(0.6), width = 0.2,
                  data = filter(combined, method == "RF\n(nonlinear)")) +
    scale_fill_manual(values = c("varpart\n(linear)" = "grey60",
                                  "RF\n(nonlinear)" = "#0E5A78")) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                       labels = scales::percent_format(accuracy = 0.1)) +
    labs(title = title,
         x = NULL, y = "Variance explained / Importance",
         fill = "Method") +
    theme_nature_pub() +
    theme(legend.position = "bottom")
}

# ── RF 变量级重要性柱状图 ──────────────────────────────────────────────
plot_rf_variable_importance <- function(rf_summary, top_n = 15,
                                         title = "RF permutation importance") {
  rf_summary |>
    mutate(variable = fct_reorder(variable, mean)) |>
    slice_max(mean, n = top_n) |>
    ggplot(aes(x = variable, y = mean)) +
    geom_col(fill = "#0E5A78", alpha = 0.85, width = 0.55) +
    geom_errorbar(aes(ymin = q_lo, ymax = q_hi), width = 0.2, linewidth = 0.3) +
    coord_flip() +
    labs(title = title, x = NULL, y = "Permutation importance (mean ± 95% CI)") +
    theme_nature_pub()
}

message("[utils_importance] loaded")
