#!/usr/bin/env Rscript
## 04b_quick_thin.R -- 极速精简：只加载 chain1，thin psi.samples 到 400 draws，gzip 保存

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tibble) })

CODE_V3 <- Sys.getenv("V3_CODE_DIR", file.path("~", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
P <- ensure_v3_dirs()

run_label <- Sys.getenv("V3_RUN_LABEL", RUN_LABEL)
psi_max <- PSI_MAX_DRAWS  # 400

chain_file <- v3_file("derived", paste0("tMsPGOcc_fit_", run_label, "_chain1"), "rds")
out_file   <- v3_file("derived", paste0("tMsPGOcc_fit_", run_label), "rds")

if (!file.exists(chain_file)) stop("Chain 1 not found: ", chain_file)
if (file.exists(out_file)) {
  message("[04b-quick] Output already exists: ", out_file)
  quit(status = 0)
}

message("[04b-quick] Loading chain 1: ", basename(chain_file))
fit <- readRDS(chain_file)

# Thin psi.samples
n_draws <- dim(fit$psi.samples)[1]
message("[04b-quick] psi.samples draws: ", n_draws)
if (n_draws > psi_max) {
  keep_idx <- round(seq(1, n_draws, length.out = psi_max))
  fit$psi.samples <- fit$psi.samples[keep_idx, , , , drop = FALSE]
  message("[04b-quick] Thinned to ", psi_max, " draws")
}

# Drop z.samples to save space
fit$z.samples <- NULL
fit$psi_max_draws <- min(n_draws, psi_max)

# Save with gzip (much faster than xz)
message("[04b-quick] Saving to: ", basename(out_file))
saveRDS(fit, out_file, compress = "gzip")
message("[04b-quick] DONE. File size: ", format(structure(file.info(out_file)$size, class = "object_size"), units = "auto"))
