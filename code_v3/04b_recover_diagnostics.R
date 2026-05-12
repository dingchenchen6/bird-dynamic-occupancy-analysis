#!/usr/bin/env Rscript
## 04b_recover_diagnostics.R  —  v3 模型诊断恢复
##
## 支持两种模式：
##   1. 单文件模式：stMsPGOcc_fit_{run_label}.rds 或 tMsPGOcc_fit_{run_label}.rds
##   2. 多链模式：stMsPGOcc_fit_{run_label}_chain{1-4}.rds 或 tMsPGOcc_fit_{run_label}_chain{1-4}.rds
## 优先查找 stMsPGOcc，回退到 tMsPGOcc
##
## 多链模式下：
##   - 加载各链独立结果
##   - 合并后验样本（沿链维度 abind）
##   - 计算 R-hat 和 ESS 诊断
##   - 保存合并后的完整 fit 对象

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_diagnostics.R"))
P <- ensure_v3_dirs()

log_time("04b", "Starting diagnostics recovery")

# ── 1. 加载模型 ──────────────────────────────────────────────────────
is_pilot <- Sys.getenv("V3_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL
n_chains_expected <- FULL_N_CHAINS

# 尝试加载：优先多链文件，回退到单文件
# 兼容 stMsPGOcc 和 tMsPGOcc 两种文件名前缀
derived_dir <- dirname(v3_file("derived", "stMsPGOcc_fit_", "rds"))
chain_pattern_st <- paste0("stMsPGOcc_fit_", run_label, "_chain[0-9]+\\.rds$")
chain_pattern_t  <- paste0("tMsPGOcc_fit_",  run_label, "_chain[0-9]+\\.rds$")

chain_files <- list.files(derived_dir, pattern = chain_pattern_st, full.names = TRUE)
if (length(chain_files) == 0) {
  # fallback: 搜索 tMsPGOcc 链文件
  chain_files <- list.files(derived_dir, pattern = chain_pattern_t, full.names = TRUE)
  if (length(chain_files) > 0) {
    message("[04b] stMsPGOcc chain files not found, using tMsPGOcc chain files")
  }
}
# 按链ID排序
chain_ids <- as.integer(gsub(".*chain([0-9]+)\\.rds$", "\\1", basename(chain_files)))
chain_files <- chain_files[order(chain_ids)]

single_file <- v3_file("derived", paste0("stMsPGOcc_fit_", run_label), "rds")
single_file_t <- v3_file("derived", paste0("tMsPGOcc_fit_", run_label), "rds")

if (length(chain_files) >= 2) {
  # ── 多链模式：加载并合并 ─────────────────────────────────────────
  message(sprintf("[04b] Found %d chain files. Combining...", length(chain_files)))

  chain_fits <- lapply(chain_files, function(f) {
    message(sprintf("[04b]   Loading: %s", basename(f)))
    readRDS(f)
  })

  # 合并后验样本
  fit <- combine_chains(chain_fits)
  message(sprintf("[04b] Combined %d chains successfully", length(chain_fits)))

  # 保存合并后的完整 fit 对象（始终保存为 stMsPGOcc 名）
  saveRDS(fit, single_file, compress = "xz")
  message(sprintf("[04b] Combined fit saved to: %s", single_file))

} else if (file.exists(single_file)) {
  # ── 单文件模式（stMsPGOcc）────────────────────────────────────────
  message(sprintf("[04b] Loading single fit file: %s", single_file))
  fit <- readRDS(single_file)
} else if (file.exists(single_file_t)) {
  # ── 单文件模式（tMsPGOcc fallback）────────────────────────────────
  message(sprintf("[04b] Loading single fit file: %s", single_file_t))
  fit <- readRDS(single_file_t)
} else {
  stop("Model fit not found (tried stMsPGOcc and tMsPGOcc). Run 04 or 04a first.")
}

# ── 2. 群落级系数摘要 ────────────────────────────────────────────────
beta_comm <- summarise_post(fit$beta.comm.samples)
alpha_comm <- summarise_post(fit$alpha.comm.samples)
tau_sq <- summarise_post(fit$tau.sq.samples)

# 空间参数
if (!is.null(fit$phi.samples)) {
  phi_summary <- summarise_post(fit$phi.samples)
  sigma_sq_summary <- summarise_post(fit$sigma.sq.samples)
  message(sprintf("[04b] Spatial range (phi): %.2f [%.2f, %.2f]",
                   phi_summary$mean, phi_summary$q025, phi_summary$q975))
  message(sprintf("[04b] Spatial variance (sigma.sq): %.4f [%.4f, %.4f]",
                   sigma_sq_summary$mean, sigma_sq_summary$q025, sigma_sq_summary$q975))
}

# ── 3. 种级系数摘要 ──────────────────────────────────────────────────
n_sp <- nrow(fit$y)
species <- rownames(fit$y) %||% paste0("sp", seq_len(n_sp))

beta_sp_list <- list()
for (i in seq_len(n_sp)) {
  beta_i <- summarise_post(fit$beta.samples[[i]])
  beta_i$species <- species[i]
  beta_sp_list[[i]] <- beta_i
}
beta_sp <- bind_rows(beta_sp_list)

# ── 4. 收敛诊断（R-hat 和 ESS）──────────────────────────────────────
diag_summary <- compute_convergence_diagnostics(fit, n_chains_expected)

# ── 5. WAIC ──────────────────────────────────────────────────────────
waic_val <- tryCatch({
  waic(fit)$waic
}, error = function(e) {
  message("[04b] WAIC computation failed: ", conditionMessage(e))
  NA_real_
})
message(sprintf("[04b] WAIC = %.2f", waic_val))

# ── 6. 保存摘要 ──────────────────────────────────────────────────────
write_csv(beta_comm |> mutate(param = "beta.comm"),
          v3_file("results", paste0("table_beta_community_", run_label)))
write_csv(alpha_comm |> mutate(param = "alpha.comm"),
          v3_file("results", paste0("table_alpha_community_", run_label)))
write_csv(tau_sq |> mutate(param = "tau.sq"),
          v3_file("results", paste0("table_tau_sq_", run_label)))
write_csv(beta_sp,
          v3_file("results", paste0("table_beta_species_", run_label)))

# 诊断摘要
if (!is.null(diag_summary)) {
  write_csv(diag_summary,
            v3_file("results", paste0("table_convergence_diagnostics_", run_label)))
  # 报告关键指标
  if ("rhat" %in% names(diag_summary)) {
    max_rhat <- max(diag_summary$rhat, na.rm = TRUE)
    min_ess  <- min(diag_summary$ess, na.rm = TRUE)
    n_above  <- sum(diag_summary$rhat > RHAT_THRESHOLD, na.rm = TRUE)
    n_below  <- sum(diag_summary$ess < ESS_THRESHOLD, na.rm = TRUE)
    message(sprintf("[04b] Convergence: max R-hat = %.3f (%d above %.2f), min ESS = %d (%d below %d)",
                    max_rhat, n_above, RHAT_THRESHOLD, min_ess, n_below, ESS_THRESHOLD))
    if (max_rhat > 1.10) {
      warning("[04b] ⚠️ Some parameters have R-hat > 1.10 — consider longer chains!")
    } else if (max_rhat > RHAT_THRESHOLD) {
      warning("[04b] ⚠️ Some parameters have R-hat > ", RHAT_THRESHOLD,
              " — may need more iterations")
    } else {
      message("[04b] ✅ All R-hat ≤ ", RHAT_THRESHOLD, " — chains have converged")
    }
  }
}

log_time("04b", sprintf("DONE: WAIC=%.2f", waic_val))

# ── 辅助函数 ──────────────────────────────────────────────────────────

#' 合并多链 stMsPGOcc 结果
#' @param chain_fits list of stMsPGOcc fit objects (each n.chains=1)
#' @return combined fit object (compatible with spOccupancy post-processing)
combine_chains <- function(chain_fits) {
  if (length(chain_fits) == 0) stop("No chain fits to combine")
  if (length(chain_fits) == 1) return(chain_fits[[1]])

  combined <- chain_fits[[1]]  # start with first chain

  # 需要沿链维度合并的后验样本 slot
  # spOccupancy 中 1-chain 结果是 matrix/vector，多链是 3D array
  sample_slots <- c(
    "beta.comm.samples",    # community occupancy coefficients
    "alpha.comm.samples",   # community detection coefficients
    "tau.sq.beta.samples",  # community variance (occ)
    "tau.sq.alpha.samples", # community variance (det)
    "tau.sq.samples",       # residual variance
    "phi.samples",          # spatial decay
    "sigma.sq.samples",     # spatial variance
    "rho.samples"           # AR1 correlation
  )

  for (slot in sample_slots) {
    if (slot %in% names(combined)) {
      combined[[slot]] <- abind_chains(
        lapply(chain_fits, function(f) f[[slot]]),
        slot
      )
    }
  }

  # Species-level slots (lists)
  sp_slots <- c("beta.samples", "alpha.samples")
  for (slot in sp_slots) {
    if (slot %in% names(combined)) {
      n_sp <- length(combined[[slot]])
      for (i in seq_len(n_sp)) {
        combined[[slot]][[i]] <- abind_chains(
          lapply(chain_fits, function(f) f[[slot]][[i]]),
          slot
        )
      }
    }
  }

  # z.samples 和 psi.samples: species × site × time × sample × chain
  # 如果存在（大对象，可能被省略）
  for (slot in c("z.samples", "psi.samples")) {
    if (slot %in% names(combined)) {
      combined[[slot]] <- abind_chains(
        lapply(chain_fits, function(f) f[[slot]]),
        slot
      )
    }
  }

  # w.samples (spatial random effects): very large, skip if absent
  # theta.samples, lambda.samples (factor model)
  for (slot in c("theta.samples", "lambda.samples")) {
    if (slot %in% names(combined)) {
      combined[[slot]] <- abind_chains(
        lapply(chain_fits, function(f) f[[slot]]),
        slot
      )
    }
  }

  # Update n.chains attribute
  combined$n.chains <- length(chain_fits)
  combined$chain_ids <- seq_len(length(chain_fits))

  combined
}

#' 沿新维度 abind 多链后验样本
#' spOccupancy single-chain: matrix [param, sample] or vector [sample]
#' Multi-chain: array [param, sample, chain] or matrix [sample, chain]
abind_chains <- function(chain_samples, slot_name = "") {
  # Remove NULL entries
  chain_samples <- chain_samples[!sapply(chain_samples, is.null)]
  if (length(chain_samples) == 0) return(NULL)
  if (length(chain_samples) == 1) {
    # Add chain dimension
    s <- chain_samples[[1]]
    if (is.vector(s)) {
      return(matrix(s, ncol = 1, dimnames = list(NULL, "chain1")))
    } else if (is.matrix(s)) {
      return(array(s, dim = c(nrow(s), ncol(s), 1),
                   dimnames = list(rownames(s), colnames(s), "chain1")))
    } else if (is.array(s) && length(dim(s)) == 3) {
      return(abind::abind(s, along = 3))  # already 3D, add chain dim
    }
    return(s)
  }

  # Check dimensions consistency
  dims <- lapply(chain_samples, function(s) {
    if (is.null(s)) return(NULL)
    dim(s)
  })
  dims <- dims[!sapply(dims, is.null)]

  if (length(dims) == 0) return(NULL)

  # All vectors → cbind into matrix [sample, chain]
  if (all(sapply(dims, is.null))) {
    return(do.call(cbind, chain_samples))
  }

  # All matrices → abind into [param, sample, chain]
  if (all(sapply(dims, length) == 2)) {
    # Check same nrow
    n_rows <- unique(sapply(dims, `[`, 1))
    if (length(n_rows) == 1) {
      result <- abind::abind(chain_samples, along = 3)
      dimnames(result)[[3]] <- paste0("chain", seq_along(chain_samples))
      return(result)
    }
  }

  # All 3D arrays → abind into [dim1, dim2, dim3, chain]
  if (all(sapply(dims, length) == 3)) {
    result <- abind::abind(chain_samples, along = 4)
    dimnames(result)[[4]] <- paste0("chain", seq_along(chain_samples))
    return(result)
  }

  # Fallback: try simple abind
  tryCatch({
    result <- abind::abind(chain_samples, along = length(dims[[1]]) + 1)
    return(result)
  }, error = function(e) {
    warning(sprintf("[04b] Could not combine slot '%s': %s", slot_name, conditionMessage(e)))
    return(chain_samples[[1]])
  })
}

#' 计算收敛诊断（R-hat 和 ESS）
compute_convergence_diagnostics <- function(fit, n_chains) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    message("[04b] coda package not available. Skipping convergence diagnostics.")
    return(NULL)
  }

  diag_list <- list()

  # Beta community samples
  if (!is.null(fit$beta.comm.samples)) {
    beta_diag <- diagnose_param_array(fit$beta.comm.samples, "beta.comm", n_chains)
    if (!is.null(beta_diag)) diag_list <- c(diag_list, list(beta_diag))
  }

  # Alpha community samples
  if (!is.null(fit$alpha.comm.samples)) {
    alpha_diag <- diagnose_param_array(fit$alpha.comm.samples, "alpha.comm", n_chains)
    if (!is.null(alpha_diag)) diag_list <- c(diag_list, list(alpha_diag))
  }

  # Phi (spatial decay)
  if (!is.null(fit$phi.samples)) {
    phi_diag <- diagnose_param_array(fit$phi.samples, "phi", n_chains)
    if (!is.null(phi_diag)) diag_list <- c(diag_list, list(phi_diag))
  }

  # Sigma.sq (spatial variance)
  if (!is.null(fit$sigma.sq.samples)) {
    sigma_diag <- diagnose_param_array(fit$sigma.sq.samples, "sigma.sq", n_chains)
    if (!is.null(sigma_diag)) diag_list <- c(diag_list, list(sigma_diag))
  }

  # Rho (AR1)
  if (!is.null(fit$rho.samples)) {
    rho_diag <- diagnose_param_array(fit$rho.samples, "rho", n_chains)
    if (!is.null(rho_diag)) diag_list <- c(diag_list, list(rho_diag))
  }

  if (length(diag_list) > 0) {
    bind_rows(diag_list)
  } else {
    NULL
  }
}

#' 诊断参数数组的 R-hat 和 ESS
diagnose_param_array <- function(samples, param_name, n_chains) {
  tryCatch({
    # 如果有链维度（最后一维），按链拆分
    if (is.array(samples) && length(dim(samples)) >= 2) {
      ndim <- length(dim(samples))

      if (ndim == 2) {
        # matrix [param, sample] — 单链
        n_params <- nrow(samples)
        results <- list()
        for (i in seq_len(n_params)) {
          chain_samples <- list(samples[i, ])
          rhat_val <- 1  # 单链无法计算 R-hat
          ess_val <- coda::effectiveSize(chain_samples[[1]])
          param_label <- if (!is.null(rownames(samples))) rownames(samples)[i] else paste0(param_name, "[", i, "]")
          results[[i]] <- tibble(
            parameter = param_label,
            group = param_name,
            rhat = rhat_val,
            ess = as.numeric(ess_val)
          )
        }
        return(bind_rows(results))
      }

      if (ndim == 3) {
        # 3D array [param, sample, chain]
        n_params <- dim(samples)[1]
        n_ch <- dim(samples)[3]
        results <- list()
        for (i in seq_len(n_params)) {
          # 提取各链样本
          chain_list <- lapply(seq_len(n_ch), function(ch) samples[i, , ch])
          # R-hat (Gelman-Rubin)
          rhat_val <- compute_rhat(chain_list)
          # ESS (total across chains)
          ess_val <- sum(sapply(chain_list, coda::effectiveSize))
          param_label <- if (!is.null(dimnames(samples)[[1]])) dimnames(samples)[[1]][i] else paste0(param_name, "[", i, "]")
          results[[i]] <- tibble(
            parameter = param_label,
            group = param_name,
            rhat = rhat_val,
            ess = as.numeric(ess_val)
          )
        }
        return(bind_rows(results))
      }

      if (ndim == 4) {
        # 4D array [dim1, dim2, dim3, chain] — e.g., species-level
        # 只诊断第一维（简化）
        n_dim1 <- dim(samples)[1]
        results <- list()
        for (i in seq_len(min(n_dim1, 20))) {  # 限制诊断数量
          chain_list <- lapply(seq_len(dim(samples)[4]), function(ch) {
            as.vector(samples[i, , , ch])
          })
          rhat_val <- compute_rhat(chain_list)
          ess_val <- sum(sapply(chain_list, coda::effectiveSize))
          results[[i]] <- tibble(
            parameter = paste0(param_name, "[", i, "]"),
            group = param_name,
            rhat = rhat_val,
            ess = as.numeric(ess_val)
          )
        }
        return(bind_rows(results))
      }
    }

    # vector: single parameter
    if (is.vector(samples)) {
      ess_val <- coda::effectiveSize(samples)
      return(tibble(
        parameter = param_name,
        group = param_name,
        rhat = 1,  # 单链
        ess = as.numeric(ess_val)
      ))
    }

    NULL
  }, error = function(e) {
    message(sprintf("[04b] Diagnostics failed for %s: %s", param_name, conditionMessage(e)))
    NULL
  })
}

#' 计算 Gelman-Rubin R-hat
compute_rhat <- function(chain_list) {
  if (length(chain_list) < 2) return(1)  # 单链无法计算

  # 确保每条链长度相同
  min_len <- min(sapply(chain_list, length))
  chains_matrix <- sapply(chain_list, function(ch) ch[1:min_len])

  n <- min_len
  m <- length(chain_list)

  # Chain means and variances
  chain_means <- colMeans(chains_matrix)
  chain_vars <- apply(chains_matrix, 2, var)

  # Between-chain variance
  B <- n * var(chain_means)
  # Within-chain variance
  W <- mean(chain_vars)

  if (W == 0) return(1)

  # Pooled variance estimate
  var_hat <- ((n - 1) / n) * W + (1 / n) * B

  # R-hat
  rhat <- sqrt(var_hat / W)
  max(rhat, 1)  # R-hat should be >= 1
}
