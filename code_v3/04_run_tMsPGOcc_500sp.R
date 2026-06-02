#!/usr/bin/env Rscript
## 04_run_tMsPGOcc_500sp.R — 500 species temporal multi-species occupancy (NO spatial)
##
## Scientific question / 科学问题:
## How do 500 bird species' occupancies change over time across China?
## 500 种鸟类的占有率如何随时间在中国变化？
##
## Strategy / 策略:
## Use tMsPGOcc (temporal, non-spatial) to avoid memory crash from spatial NNGP + 15 factors.
## 使用 tMsPGOcc（仅时间，无空间）以避免 spatial NNGP + 15 factors 导致的内存崩溃。
##
## Workflow / 流程:
## 1. Load survey + environment data
## 2. Filter to top 500 species
## 3. Fit tMsPGOcc with AR1 + 10 factors
## 4. Save per-chain RDS

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(spOccupancy)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
P <- ensure_v3_dirs()

# Chain ID from env var
chain_id <- as.integer(Sys.getenv("V3_CHAIN_ID", NA))
n_omp <- as.integer(Sys.getenv("V3_N_OMP", "32"))

# Force 500sp / tMsPGOcc settings
RUN_LABEL_T <- "v3_full_500sp_ar1_temporal"
MAX_N_SP <- 500L
N_FACTORS_T <- 10L  # Reduced from 15 (no spatial = simpler model)
N_BATCH <- 600L
N_BURN  <- 8000L
N_THIN  <- 2L
N_CHAINS <- 1L  # 1 chain per process (orchestrated externally)

log_time("04_t500", sprintf("Starting tMsPGOcc 500sp (chain %s, n_omp=%d)",
                            ifelse(is.na(chain_id), "single", chain_id), n_omp))

# ── 1. Load data ──────────────────────────────────────────────────────
survey <- safe_read(v3_file("derived", paste0("survey_history", GRID_TAG, "_v3"), "rds"))
if (is.null(survey)) {
  survey <- safe_read(v3_file("derived", "survey_history_v3", "rds"))
}
grid_env <- safe_read(v3_file("derived", paste0("grid_environment", GRID_TAG, "_v3"), "rds"))
if (is.null(grid_env)) {
  grid_env <- safe_read(v3_file("derived", "grid_environment_v3", "rds"))
}

# ── 2. Pick top 500 species by detections ─────────────────────────────
candidates_all <- read_csv_safe(v3_file("results", "table_candidate_species_v3.csv"))
if (!is.null(candidates_all)) {
  candidate_species <- head(candidates_all$species, MAX_N_SP)
} else {
  candidate_species <- head(survey$species, MAX_N_SP)
}
n_sp <- length(candidate_species)
message(sprintf("[04_t500] Selected %d species", n_sp))

# ── 3. Build y array: species x site x period x rep ───────────────────
y_full <- survey$y
sp_idx <- which(survey$species %in% candidate_species)
y_sub <- y_full[sp_idx, , , , drop = FALSE]
dimnames(y_sub)[[1]] <- survey$species[sp_idx]
n_sp <- dim(y_sub)[1]; n_sites <- dim(y_sub)[2]
n_periods <- dim(y_sub)[3]; n_rep <- dim(y_sub)[4]
message(sprintf("[04_t500] y dims: %d sp x %d sites x %d periods x %d reps",
                n_sp, n_sites, n_periods, n_rep))

# ── 4. Detection covariates ───────────────────────────────────────────
det_cov <- safe_read(v3_file("derived", "detection_covariates_v3", "rds"))
if (is.null(det_cov)) stop("Detection covariates not found.")

# ── 5. Occupancy covariates from grid_env ────────────────────────────
grid_env <- grid_env[match(survey$sites, grid_env$grid_cell), ]
occ_vars <- c("bio4", "bio7", "bio11", "bio13", "elev_mean", "elev_sd",
              "texture_shannon", "habitat_diversity_shannon", "hfi_mean",
              "landcover_built", "landcover_cropland",
              "centroid_lon", "centroid_lat")
occ_cov_df <- as.data.frame(grid_env[, occ_vars])
# Fill NAs with column median
for (cc in occ_vars) {
  v <- occ_cov_df[[cc]]
  if (any(is.na(v))) v[is.na(v)] <- median(v, na.rm = TRUE)
  occ_cov_df[[cc]] <- as.numeric(scale(v))
}
# Add year_scaled (time covariate, replicated per period)
year_scaled_per_period <- scale(seq_len(n_periods))[, 1]
occ_cov_list <- list()
for (v in occ_vars) {
  occ_cov_list[[v]] <- matrix(rep(occ_cov_df[[v]], n_periods),
                              nrow = n_sites, ncol = n_periods)
}
occ_cov_list[["year_scaled"]] <- matrix(rep(year_scaled_per_period, each = n_sites),
                                         nrow = n_sites, ncol = n_periods)

# ── 6. Build data_list for tMsPGOcc ──────────────────────────────────
data_list <- list(
  y = y_sub,
  occ.covs = occ_cov_list,
  det.covs = det_cov
)

# Drop sites with no detections
site_has_data <- apply(y_sub, 2, function(s) any(!is.na(s)))
if (sum(!site_has_data) > 0) {
  message(sprintf("[04_t500] Dropping %d sites with no detections", sum(!site_has_data)))
  data_list$y <- data_list$y[, site_has_data, , , drop = FALSE]
  for (v in names(data_list$occ.covs)) {
    data_list$occ.covs[[v]] <- data_list$occ.covs[[v]][site_has_data, , drop = FALSE]
  }
  for (v in names(data_list$det.covs)) {
    data_list$det.covs[[v]] <- data_list$det.covs[[v]][site_has_data, , , drop = FALSE]
  }
}

n_sites_use <- dim(data_list$y)[2]
message(sprintf("[04_t500] Final dims: %d sp x %d sites x %d periods", n_sp, n_sites_use, n_periods))

# ── 7. Priors / inits / tuning ────────────────────────────────────────
occ_formula <- as.formula(OCC_FORMULA_STR)
det_formula <- as.formula(DET_FORMULA_STR)

priors <- list(
  beta.comm.normal = list(mean = 0, var = 2.72),
  alpha.comm.normal = list(mean = 0, var = 2.72),
  tau.sq.beta.ig = list(a = 0.1, b = 0.1),
  tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
  sigma.sq.t.ig = list(a = 2, b = 0.5),
  rho.unif = list(a = -1, b = 1)
)

inits <- list(
  beta.comm = 0, alpha.comm = 0,
  tau.sq.beta = 1, tau.sq.alpha = 1,
  sigma.sq.t = 0.5, rho = 0
)

# ── 8. Run tMsPGOcc ───────────────────────────────────────────────────
message(sprintf("[04_t500] Running tMsPGOcc: %s | n.factors=%d | n.batch=%d | n.burn=%d | n.omp=%d",
                RUN_LABEL_T, N_FACTORS_T, N_BATCH, N_BURN, n_omp))

t_start <- Sys.time()

fit <- tMsPGOcc(
  occ.formula   = occ_formula,
  det.formula   = det_formula,
  data          = data_list,
  inits         = inits,
  priors        = priors,
  n.factors     = N_FACTORS_T,
  n.batch       = N_BATCH,
  batch.length  = 25,
  accept.rate   = 0.43,
  n.omp.threads = n_omp,
  verbose       = TRUE,
  ar1           = TRUE,
  n.report      = max(1, floor(N_BATCH / 10)),
  n.burn        = N_BURN,
  n.thin        = N_THIN,
  n.chains      = N_CHAINS
)

t_elapsed <- difftime(Sys.time(), t_start, units = "mins")
message(sprintf("[04_t500] Model fit completed in %.1f minutes", as.numeric(t_elapsed)))

# ── 9. Save ───────────────────────────────────────────────────────────
if (!is.na(chain_id)) {
  fit_path <- v3_file("derived", paste0("tMsPGOcc_fit_", RUN_LABEL_T, "_chain", chain_id), "rds")
} else {
  fit_path <- v3_file("derived", paste0("tMsPGOcc_fit_", RUN_LABEL_T), "rds")
}

gc()
saveRDS(fit, fit_path)
message(sprintf("[04_t500] Saved: %s (%.1f MB)", fit_path, file.size(fit_path) / 1024^2))

log_time("04_t500", "DONE")
