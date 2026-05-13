#!/usr/bin/env Rscript
## utils_diagnostics.R  —  v3 诊断工具模块
##
## 提供种级收敛诊断、WAIC 计算、DHARMa gate 检查

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

# ── 种级收敛诊断提取 ─────────────────────────────────────────────────
#' 从 stMsPGOcc / tMsPGOcc 输出中提取物种水平 R-hat 和 ESS
#'
#' @param fit 模型拟合对象（spOccupancy 输出）
#' @param n_species 物种数
#' @param species_names 物种名向量
#' @return tibble: species, param, rhat, ess
extract_species_diag <- function(fit, n_species = NULL, species_names = NULL) {
  # 尝试获取 Rhat / ESS
  rhat_out <- fit$rhat
  ess_out  <- fit$ESS

  if (is.null(rhat_out) || is.null(ess_out)) {
    warning("Model fit missing rhat/ESS slots. Try running summary().")
    return(NULL)
  }

  # 确定物种数
  if (is.null(n_species)) {
    n_species <- if (!is.null(fit$y)) nrow(fit$y) else NROW(rhat_out)
  }
  if (is.null(species_names)) {
    species_names <- if (!is.null(rownames(fit$y))) {
      rownames(fit$y)
    } else {
      paste0("sp", seq_len(n_species))
    }
  }

  # 构建诊断表
  diag_list <- list()
  param_names <- rownames(rhat_out) %||% paste0("param", seq_len(nrow(rhat_out)))

  for (p in seq_along(param_names)) {
    diag_list[[p]] <- tibble(
      species = species_names,
      param   = param_names[p],
      rhat    = as.numeric(rhat_out[p, ]),
      ess     = as.numeric(ess_out[p, ])
    )
  }

  bind_rows(diag_list)
}

# ── WAIC 计算 ────────────────────────────────────────────────────────
#' 从 spOccupancy 模型计算 WAIC
#'
#' @param fit 模型拟合对象
#' @return 数值向量: waic, p_waic
compute_waic <- function(fit) {
  if (!requireNamespace("spOccupancy", quietly = TRUE)) {
    stop("spOccupancy package required for WAIC computation")
  }

  # spOccupancy 自带 WAIC
  if ("waic" %in% names(fit)) {
    return(fit$waic)
  }

  # 若无预计算，使用 spOccupancy::waicOcc
  if (requireNamespace("spOccupancy", quietly = TRUE) &&
      exists("waicOcc", envir = asNamespace("spOccupancy"))) {
    return(spOccupancy::waicOcc(fit))
  }

  warning("WAIC not available in fit object. Run waicOcc() manually.")
  NA_real_
}

# ── DHARMa gate 检查 ──────────────────────────────────────────────────
#' 对 brms 模型执行 DHARMa 残差诊断，自动判断是否通过
#'
#' @param brms_fit brms 模型对象
#' @param save_dir 保存目录
#' @param stem 文件名前缀
#' @return tibble: test, statistic, p_value, passed
check_dharma_gate <- function(brms_fit, save_dir = NULL, stem = "dharma") {
  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    warning("DHARMa package not installed. Skipping gate check.")
    return(NULL)
  }

  # 模拟残差（brms 无随机效应时 coef() 可能报错，用 tryCatch 保护）
  sim_res <- tryCatch(
    DHARMa::simulateResiduals(brms_fit, n = 250, plot = FALSE),
    error = function(e) {
      warning("[dharma_gate] simulateResiduals failed: ", e$message, call. = FALSE)
      NULL
    }
  )
  if (is.null(sim_res)) return(NULL)

  # 过离散检验
  od_test <- DHARMa::testDispersion(sim_res, plot = FALSE)
  # 均匀性检验
  unif_test <- DHARMa::testUniformity(sim_res, plot = FALSE)

  # 异常值检验
  out_test <- DHARMa::testOutliers(sim_res, plot = FALSE)

  result <- tibble(
    test     = c("overdispersion", "uniformity", "outliers"),
    statistic = c(od_test$statistic, unif_test$statistic, out_test$statistic),
    p_value  = c(od_test$p.value, unif_test$p.value, out_test$p.value),
    passed   = c(od_test$p.value > 0.05, unif_test$p.value > 0.05,
                 out_test$p.value > 0.05)
  )

  # 保存诊断图
  if (!is.null(save_dir)) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    png_path <- file.path(save_dir, paste0(stem, ".png"))
    png(png_path, width = 10, height = 5, units = "in", res = 300)
    par(mfrow = c(1, 2))
    DHARMa::plotSimulatedResiduals(sim_res)
    dev.off()
    message("[dharma] saved → ", png_path)
  }

  # Gate 结果
  all_passed <- all(result$passed, na.rm = TRUE)
  if (!all_passed) {
    failed <- result$test[!result$passed]
    warning("[dharma_gate] FAILED: ", paste(failed, collapse = ", "),
            call. = FALSE)
  } else {
    message("[dharma_gate] PASSED all tests")
  }

  result
}

# ── MCMC 收敛汇总 ─────────────────────────────────────────────────────
#' 汇总种级和群落级收敛诊断
summarise_convergence <- function(diag_df, rhat_thresh = 1.05,
                                   ess_thresh = 200) {
  diag_df |>
    group_by(param) |>
    summarise(
      n_species       = n(),
      rhat_max        = max(rhat, na.rm = TRUE),
      rhat_gt_thresh  = sum(rhat > rhat_thresh, na.rm = TRUE),
      ess_min         = min(ess, na.rm = TRUE),
      ess_lt_thresh   = sum(ess < ess_thresh, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      rhat_ok = rhat_gt_thresh == 0,
      ess_ok  = ess_lt_thresh == 0,
      overall = rhat_ok & ess_ok
    )
}

message("[utils_diagnostics] loaded")
