#!/usr/bin/env Rscript
## 04b_recover_diagnostics_light.R -- 轻量版诊断恢复
## 策略：加载每链后仅保留 psi.samples 的 400 draws 子集，不保存 z.samples
## 合并后 fit 对象从 ~60GB 降到 ~5GB，加载/保存快 10 倍

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_diagnostics.R"))
P <- ensure_v3_dirs()

log_time("04b", "Starting LIGHT diagnostics recovery")

is_pilot <- Sys.getenv("V3_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL
n_chains_expected <- FULL_N_CHAINS
psi_max <- PSI_MAX_DRAWS  # 400

# ============================================================
# 辅助函数
# ============================================================

abind_chains <- function(chain_samples, slot_name = "") {
  chain_samples <- chain_samples[!sapply(chain_samples, is.null)]
  if (length(chain_samples) == 0) return(NULL)
  if (length(chain_samples) == 1) {
    s <- chain_samples[[1]]
    if (is.vector(s)) {
      return(matrix(s, ncol = 1, dimnames = list(NULL, "chain1")))
    } else if (is.matrix(s)) {
      return(array(s, dim = c(nrow(s), ncol(s), 1),
                   dimnames = list(rownames(s), colnames(s), "chain1")))
    } else if (is.array(s) && length(dim(s)) == 3) {
      return(s)
    }
    return(s)
  }
  dims <- lapply(chain_samples, function(s) { if (is.null(s)) return(NULL); dim(s) })
  dims <- dims[!sapply(dims, is.null)]
  if (length(dims) == 0) return(NULL)
  if (all(sapply(dims, is.null))) return(do.call(cbind, chain_samples))
  if (all(sapply(dims, length) == 2)) {
    n_rows <- unique(sapply(dims, `[`, 1))
    if (length(n_rows) == 1) {
      result <- abind::abind(chain_samples, along = 3)
      dimnames(result)[[3]] <- paste0("chain", seq_along(chain_samples))
      return(result)
    }
  }
  if (all(sapply(dims, length) == 3)) {
    result <- abind::abind(chain_samples, along = 4)
    dimnames(result)[[4]] <- paste0("chain", seq_along(chain_samples))
    return(result)
  }
  tryCatch({
    result <- abind::abind(chain_samples, along = length(dims[[1]]) + 1)
    return(result)
  }, error = function(e) {
    warning(sprintf("[04b] Could not combine slot '%s': %s", slot_name, conditionMessage(e)))
    return(chain_samples[[1]])
  })
}

compute_rhat <- function(chain_list) {
  if (length(chain_list) < 2) return(1)
  min_len <- min(sapply(chain_list, length))
  chains_matrix <- sapply(chain_list, function(ch) ch[1:min_len])
  n <- min_len; m <- length(chain_list)
  chain_means <- colMeans(chains_matrix)
  chain_vars <- apply(chains_matrix, 2, var)
  B <- n * var(chain_means); W <- mean(chain_vars)
  if (W == 0) return(1)
  var_hat <- ((n - 1) / n) * W + (1 / n) * B
  max(sqrt(var_hat / W), 1)
}

diagnose_param_array <- function(samples, param_name, n_chains) {
  tryCatch({
    if (is.array(samples) && length(dim(samples)) >= 2) {
      ndim <- length(dim(samples))
      if (ndim == 2) {
        n_params <- nrow(samples)
        results <- list()
        for (i in seq_len(n_params)) {
          ess_val <- coda::effectiveSize(samples[i, ])
          param_label <- if (!is.null(rownames(samples))) rownames(samples)[i] else paste0(param_name, "[", i, "]")
          results[[i]] <- tibble(parameter = param_label, group = param_name, rhat = 1, ess = as.numeric(ess_val))
        }
        return(bind_rows(results))
      }
      if (ndim == 3) {
        n_params <- dim(samples)[1]; n_ch <- dim(samples)[3]
        results <- list()
        for (i in seq_len(n_params)) {
          chain_list <- lapply(seq_len(n_ch), function(ch) samples[i, , ch])
          rhat_val <- compute_rhat(chain_list)
          ess_val <- sum(sapply(chain_list, coda::effectiveSize))
          param_label <- if (!is.null(dimnames(samples)[[1]])) dimnames(samples)[[1]][i] else paste0(param_name, "[", i, "]")
          results[[i]] <- tibble(parameter = param_label, group = param_name, rhat = rhat_val, ess = as.numeric(ess_val))
        }
        return(bind_rows(results))
      }
    }
    if (is.vector(samples)) {
      ess_val <- coda::effectiveSize(samples)
      return(tibble(parameter = param_name, group = param_name, rhat = 1, ess = as.numeric(ess_val)))
    }
    NULL
  }, error = function(e) {
    message(sprintf("[04b] Diagnostics failed for %s: %s", param_name, conditionMessage(e)))
    NULL
  })
}

compute_convergence_diagnostics <- function(fit, n_chains) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    message("[04b] coda package not available. Skipping convergence diagnostics.")
    return(NULL)
  }
  diag_list <- list()
  for (slot_name in c("beta.comm.samples", "alpha.comm.samples", "phi.samples", "sigma.sq.samples", "rho.samples")) {
    if (!is.null(fit[[slot_name]])) {
      short_name <- gsub("\\.samples$", "", slot_name)
      d <- diagnose_param_array(fit[[slot_name]], short_name, n_chains)
      if (!is.null(d)) diag_list <- c(diag_list, list(d))
    }
  }
  if (length(diag_list) > 0) bind_rows(diag_list) else NULL
}

# ============================================================
# 主逻辑：轻量合并
# ============================================================

chain_files <- list.files(
  dirname(v3_file("derived", "tMsPGOcc_fit_", "rds")),
  pattern = paste0("tMsPGOcc_fit_", run_label, "_chain[0-9]+\\.rds$"),
  full.names = TRUE
)
chain_ids <- as.integer(gsub(".*chain([0-9]+)\\.rds$", "\\1", basename(chain_files)))
chain_files <- chain_files[order(chain_ids)]
single_file <- v3_file("derived", paste0("tMsPGOcc_fit_", run_label), "rds")

if (length(chain_files) == 0) {
  if (file.exists(single_file)) {
    message(sprintf("[04b] Loading single fit file: %s", single_file))
    fit <- readRDS(single_file)
  } else {
    stop("Model fit not found. Run 04 or 04a first.")
  }
} else {
  message(sprintf("[04b] Found %d chain files. Light merge...", length(chain_files)))

  # 加载第一个链作为基础
  message(sprintf("[04b] Loading chain 1: %s", basename(chain_files[1])))
  fit <- readRDS(chain_files[1])
  n_draws_per_chain <- dim(fit$psi.samples)[1]
  message(sprintf("[04b]   Draws per chain: %d", n_draws_per_chain))

  # 计算所有链的总draws和抽取索引
  total_draws <- n_draws_per_chain * length(chain_files)
  n_keep <- min(psi_max, total_draws)
  keep_idx <- round(seq(1, total_draws, length.out = n_keep))
  message(sprintf("[04b]   Total draws: %d, keeping: %d", total_draws, n_keep))

  # 确定每个保留索引来自哪个链和哪个draw
  chain_per_idx <- ((keep_idx - 1) %/% n_draws_per_chain) + 1
  draw_per_idx <- ((keep_idx - 1) %% n_draws_per_chain) + 1

  # 从第一个链提取需要的 draws
  keep_from_chain1 <- draw_per_idx[chain_per_idx == 1]
  psi_list <- list(fit$psi.samples[keep_from_chain1, , , , drop = FALSE])

  # 释放第一个链的大数组（但保留基础结构）
  fit$psi.samples <- NULL
  fit$z.samples <- NULL
  gc()

  # 加载其余链并提取对应 draws
  if (length(chain_files) > 1) {
    for (ch in 2:length(chain_files)) {
      message(sprintf("[04b] Loading chain %d: %s", ch, basename(chain_files[ch])))
      fit_ch <- readRDS(chain_files[ch])
      keep_from_ch <- draw_per_idx[chain_per_idx == ch]
      psi_list[[ch]] <- fit_ch$psi.samples[keep_from_ch, , , , drop = FALSE]
      rm(fit_ch); gc()
    }
  }

  # 沿draw维度合并 psi.samples
  message("[04b] Merging psi.samples...")
  fit$psi.samples <- abind::abind(psi_list, along = 1)
  dimnames(fit$psi.samples)[[1]] <- paste0("draw", seq_len(n_keep))
  rm(psi_list); gc()

  # 合并其他参数（使用 abind_chains，沿 chain 维度）
  sample_slots <- c("beta.comm.samples", "alpha.comm.samples",
                    "tau.sq.beta.samples", "tau.sq.alpha.samples",
                    "tau.sq.samples", "phi.samples", "sigma.sq.samples", "rho.samples")

  # 重新加载所有链的小参数
  message("[04b] Merging parameter samples...")
  chain_fits <- lapply(chain_files, function(f) {
    fc <- readRDS(f)
    # 只保留需要的参数，删除大数组以节省内存
    fc$psi.samples <- NULL
    fc$z.samples <- NULL
    fc
  })

  for (slot in sample_slots) {
    if (slot %in% names(fit)) {
      merged <- abind_chains(
        lapply(chain_fits, function(f) f[[slot]]), slot
      )
      if (!is.null(merged)) {
        fit[[slot]] <- merged
      } else {
        message(sprintf("[04b] WARNING: abind_chains returned NULL for %s, keeping chain1 only", slot))
      }
    }
  }

  # 合并物种级参数（处理某些链中物种缺失的情况）
  sp_slots <- c("beta.samples", "alpha.samples")
  for (slot in sp_slots) {
    if (slot %in% names(fit)) {
      n_sp <- length(fit[[slot]])
      for (i in seq_len(n_sp)) {
        samples_list <- lapply(chain_fits, function(f) {
          if (i <= length(f[[slot]]) && !is.null(f[[slot]][[i]])) f[[slot]][[i]] else NULL
        })
        samples_list <- samples_list[!sapply(samples_list, is.null)]
        if (length(samples_list) > 0) {
          merged <- abind_chains(samples_list, slot)
          if (!is.null(merged)) {
            fit[[slot]][[i]] <- merged
          } else {
            message(sprintf("[04b] WARNING: abind_chains returned NULL for %s species %d, skipping merge", slot, i))
          }
        } else {
          message(sprintf("[04b] WARNING: %s for species %d missing in all chains, keeping chain1 only", slot, i))
        }
      }
    }
  }

  rm(chain_fits); gc()

  fit$n.chains <- length(chain_files)
  fit$chain_ids <- seq_len(length(chain_files))
  fit$psi_max_draws <- n_keep
  message(sprintf("[04b] Light merge complete. psi.samples dims: %s", paste(dim(fit$psi.samples), collapse="x")))

  # 保存轻量合并文件（不使用 xz，用 gzip 更快）
  message(sprintf("[04b] Saving light fit to: %s", single_file))
  saveRDS(fit, single_file, compress = "gzip")
  message("[04b] Saved.")
}

# 参数汇总
message("[04b] Summarising community parameters...")
beta_comm <- summarise_post(fit$beta.comm.samples)
alpha_comm <- summarise_post(fit$alpha.comm.samples)
tau_sq <- summarise_post(fit$tau.sq.samples)

if (!is.null(fit$phi.samples)) {
  phi_summary <- summarise_post(fit$phi.samples)
  sigma_sq_summary <- summarise_post(fit$sigma.sq.samples)
  message(sprintf("[04b] Spatial range (phi): %.2f [%.2f, %.2f]",
                   phi_summary$mean, phi_summary$q025, phi_summary$q975))
  message(sprintf("[04b] Spatial variance (sigma.sq): %.4f [%.4f, %.4f]",
                   sigma_sq_summary$mean, sigma_sq_summary$q025, sigma_sq_summary$q975))
}

# 物种级参数
message("[04b] Summarising species parameters...")
n_sp <- nrow(fit$y)
species <- rownames(fit$y) %||% paste0("sp", seq_len(n_sp))
beta_sp_list <- list()
for (i in seq_len(n_sp)) {
  beta_i <- summarise_post(fit$beta.samples[[i]])
  beta_i$species <- species[i]
  beta_sp_list[[i]] <- beta_i
}
beta_sp <- bind_rows(beta_sp_list)

# 收敛诊断
message("[04b] Computing convergence diagnostics...")
diag_summary <- compute_convergence_diagnostics(fit, n_chains_expected)

# WAIC（轻量版 fit 可能无法计算 WAIC，因为缺少 z.samples）
message("[04b] Computing WAIC...")
waic_val <- tryCatch({
  waic(fit)$waic
}, error = function(e) {
  message("[04b] WAIC computation failed (expected without z.samples): ", conditionMessage(e))
  NA_real_
})
if (!is.na(waic_val)) message(sprintf("[04b] WAIC = %.2f", waic_val))

# 保存结果表
write_csv(beta_comm |> mutate(param = "beta.comm"),
          v3_file("results", paste0("table_beta_community_", run_label)))
write_csv(alpha_comm |> mutate(param = "alpha.comm"),
          v3_file("results", paste0("table_alpha_community_", run_label)))
write_csv(tau_sq |> mutate(param = "tau.sq"),
          v3_file("results", paste0("table_tau_sq_", run_label)))
write_csv(beta_sp,
          v3_file("results", paste0("table_beta_species_", run_label)))

if (!is.null(diag_summary)) {
  write_csv(diag_summary,
            v3_file("results", paste0("table_convergence_diagnostics_", run_label)))
  if ("rhat" %in% names(diag_summary)) {
    max_rhat <- max(diag_summary$rhat, na.rm = TRUE)
    min_ess  <- min(diag_summary$ess, na.rm = TRUE)
    n_above  <- sum(diag_summary$rhat > RHAT_THRESHOLD, na.rm = TRUE)
    n_below  <- sum(diag_summary$ess < ESS_THRESHOLD, na.rm = TRUE)
    message(sprintf("[04b] Convergence: max R-hat = %.3f (%d above %.2f), min ESS = %d (%d below %d)",
                    max_rhat, n_above, RHAT_THRESHOLD, min_ess, n_below, ESS_THRESHOLD))
    if (max_rhat > 1.10) {
      warning("[04b] Some parameters have R-hat > 1.10 -- consider longer chains!")
    } else if (max_rhat > RHAT_THRESHOLD) {
      warning("[04b] Some parameters have R-hat > ", RHAT_THRESHOLD, " -- may need more iterations")
    } else {
      message("[04b] All R-hat <= ", RHAT_THRESHOLD, " -- chains have converged")
    }
  }
}

log_time("04b", sprintf("DONE: light merge, WAIC=%s", ifelse(is.na(waic_val), "N/A", sprintf("%.2f", waic_val))))
