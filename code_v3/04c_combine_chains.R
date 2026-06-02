#!/usr/bin/env Rscript
## 04c_combine_chains.R — 4链合并 + 后半段100 draws
##
## 科学问题 / Scientific question:
##   4链合并后验是否比单链后验CI更可靠？
##   Is 4-chain merged posterior more reliable than single-chain CI?
##
## 分析目标 / Objective:
##   1. 加载4个链文件
##   2. 对psi.samples: 从每链后半段(1001-4000)均匀取100 draws
##   3. 合并为400 draws
##   4. 对其他参数: 4链abind合并(用于R-hat/ESS诊断)
##   5. 计算收敛诊断
##
## 输入 / Input:
##   data/derived_v3/tMsPGOcc_fit_*_04c_chain{1-4}.rds
##
## 输出 / Output:
##   data/derived_v3/tMsPGOcc_fit_v3_full_200sp_ar1_spatial_04c.rds
##   results_v3/table_convergence_04c.csv
##
## 主包 / Main packages: abind, readr, dplyr

suppressPackageStartupMessages({
  library(abind); library(readr); library(dplyr); library(tidyr); library(tibble)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
P <- ensure_v3_dirs()

run_label <- Sys.getenv("V3_RUN_LABEL", "v3_full_200sp_ar1_spatial_04c")
derived_dir <- DIRS$derived
n_chains <- 4L

# 后半段采样参数 / Latter-half sampling parameters
DRAWS_PER_CHAIN <- 100L
PSI_SKIP_INITIAL <- 1000L  # 跳过burn后前1000步 / skip first 1000 post-burn

message(sprintf("[04c-combine] Combining %d chains for: %s", n_chains, run_label))

# ── 1. 加载链文件 ──────────────────────────────────────────────────
chain_fits <- list()
for (ch in seq_len(n_chains)) {
  chain_file <- file.path(derived_dir,
    paste0("tMsPGOcc_fit_", run_label, "_chain", ch, ".rds"))
  if (!file.exists(chain_file)) {
    stop("Chain ", ch, " not found: ", chain_file)
  }
  fsize <- file.info(chain_file)$size
  message(sprintf("[04c-combine] Loading chain %d: %.1f GB", ch, fsize / 1e9))
  chain_fits[[ch]] <- readRDS(chain_file)
}

# ── 2. psi.samples: 4链各取后半段100 draws ──────────────────────────
n_draws_per_chain <- dim(chain_fits[[1]]$psi.samples)[1]
message(sprintf("[04c-combine] Draws per chain: %d", n_draws_per_chain))

# 计算后半段采样范围 / Calculate latter-half sampling range
if (n_draws_per_chain > PSI_SKIP_INITIAL) {
  start_idx <- PSI_SKIP_INITIAL + 1
  end_idx <- n_draws_per_chain
  draw_idx_per_chain <- round(seq(start_idx, end_idx, length.out = DRAWS_PER_CHAIN))
  message(sprintf("[04c-combine] Sampling from latter half: %d-%d → %d draws",
                  start_idx, end_idx, DRAWS_PER_CHAIN))
} else {
  # 如果链太短, 均匀取全部 / If chain too short, sample uniformly
  draw_idx_per_chain <- round(seq(1, n_draws_per_chain, length.out = DRAWS_PER_CHAIN))
  message(sprintf("[04c-combine] Chain shorter than skip, sampling uniformly: %d draws", DRAWS_PER_CHAIN))
}

psi_list <- list()
for (ch in seq_len(n_chains)) {
  psi_list[[ch]] <- chain_fits[[ch]]$psi.samples[draw_idx_per_chain, , , , drop = FALSE]
  # 释放链文件的psi以节省内存 / Free chain psi to save memory
  chain_fits[[ch]]$psi.samples <- NULL
  chain_fits[[ch]]$z.samples <- NULL
  gc()
}

# 沿draw维度合并 / Merge along draw dimension
message("[04c-combine] Merging psi.samples...")
psi_combined <- abind(psi_list, along = 1)
dimnames(psi_combined)[[1]] <- paste0("draw", seq_len(n_chains * DRAWS_PER_CHAIN))
message(sprintf("[04c-combine] psi.samples dims: %s", paste(dim(psi_combined), collapse = " × ")))
rm(psi_list); gc()

# ── 3. 其他参数: 4链按draw维度合并(用于诊断) ───────────────────────
combined <- chain_fits[[1]]
merge_draw_samples <- function(slot) {
  if (!slot %in% names(combined)) return(NULL)
  objs <- lapply(chain_fits, function(f) f[[slot]])
  if (is.null(objs[[1]])) return(NULL)
  if (is.list(objs[[1]]) && !is.data.frame(objs[[1]])) {
    n_slot <- length(objs[[1]])
    out <- objs[[1]]
    for (i in seq_len(n_slot)) {
      out[[i]] <- do.call(rbind, lapply(objs, function(o) o[[i]]))
    }
    return(out)
  }
  do.call(rbind, objs)
}

sample_slots <- c("beta.comm.samples", "alpha.comm.samples",
                  "tau.sq.beta.samples", "tau.sq.alpha.samples",
                  "tau.sq.samples", "phi.samples", "sigma.sq.samples", "rho.samples",
                  "beta.samples", "alpha.samples", "theta.samples", "lambda.samples")
for (slot in sample_slots) {
  merged_slot <- merge_draw_samples(slot)
  if (!is.null(merged_slot)) combined[[slot]] <- merged_slot
}

combined$psi.samples <- psi_combined
combined$n.chains <- n_chains
combined$chain_ids <- seq_len(n_chains)
combined$psi_merge_method <- "4chain_latter_half_100"
combined$psi_skip_initial <- PSI_SKIP_INITIAL

rm(chain_fits, psi_combined); gc()

# ── 4. 收敛诊断 ──────────────────────────────────────────────────────
source(file.path(CODE_V3, "utils_diagnostics.R"))

diag_summary <- NULL
if (exists("compute_convergence_diagnostics", mode = "function")) {
  diag_summary <- tryCatch(
    compute_convergence_diagnostics(combined, n_chains),
    error = function(e) {
      message("[04c-combine] Convergence diagnostics skipped: ", conditionMessage(e))
      NULL
    }
  )
} else {
  message("[04c-combine] Convergence diagnostics function not found; skipping.")
}

# ── 5. 保存合并文件 ──────────────────────────────────────────────────
out_file <- v3_file("derived", paste0("tMsPGOcc_fit_", run_label), "rds")
message(sprintf("[04c-combine] Saving to: %s", basename(out_file)))
saveRDS(combined, out_file, compress = "gzip")
out_size <- file.info(out_file)$size
message(sprintf("[04c-combine] Saved: %.1f GB", out_size / 1e9))

# 保存收敛诊断表 / Save convergence diagnostics
if (!is.null(diag_summary) && length(diag_summary) > 0) {
  diag_df <- tryCatch({
    bind_rows(lapply(names(diag_summary), function(nm) {
      vals <- diag_summary[[nm]]
      if (is.numeric(vals) && length(vals) > 0) {
        tibble(parameter = nm, n = length(vals),
               min_rhat = min(vals, na.rm = TRUE), max_rhat = max(vals, na.rm = TRUE),
               median_rhat = median(vals, na.rm = TRUE),
               n_above_1.05 = sum(vals > 1.05, na.rm = TRUE))
      } else NULL
    }))
  }, error = function(e) NULL)

  if (!is.null(diag_df) && nrow(diag_df) > 0) {
    write_csv(diag_df, v3_file("results", paste0("table_convergence_", run_label)))
    message("[04c-combine] Convergence diagnostics saved")
  }
}

rm(combined); gc()
message("[04c-combine] DONE")
