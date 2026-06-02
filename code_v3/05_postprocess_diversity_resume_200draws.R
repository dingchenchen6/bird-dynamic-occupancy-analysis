#!/usr/bin/env Rscript
## 05_postprocess_diversity_resume.R — Resume from brms + psi save + Section 10
## 05_postprocess_diversity.R 补完脚本：避免重跑 diversity/temporal/Baselga/MK（~17h）
##
## 科学问题 / Scientific question:
## 当 05_postprocess_diversity.R 在前半部分成功、后半部分（brms/psi/物种趋势）失败时，
## 从已保存的 CSV 结果继续执行，无需重算社区指标。
## When the first half of 05 succeeds but the second half (brms/psi/species trends) fails,
## resume from saved CSV outputs without recomputing community metrics.
##
## 输入数据 / Input:
##   - fit (tMsPGOcc/stMsPGOcc model object)
##   - grid_env
##   - table_trend_summary_*.csv (saved by previous run)
##   - trait_ext, phylo
##
## 主要流程 / Workflow:
## 1. Load config & model / 加载配置与模型
## 2. Recover trend_summary_df from CSV / 从CSV恢复趋势汇总
## 3. brms driver regression / 驱动因子回归
## 4. Save thinned psi samples / 保存 psi 抽取
## 5. Species-level trend analysis / 物种水平趋势分析
##
## 预期输出 / Expected output:
##   - brms_driver_trend_*.rds, brms_loo_trend_*.rds
##   - psi_samples_thinned_*.rds
##   - table_species_trend_*.csv, table_species_trend_classify_*.csv
##   - table_species_trend_traits_*.csv, table_trait_trend_correlation_*.csv
##   - table_species_env_trend_*.csv, table_env_trend_correlation_*.csv
##   - table_species_hotspot_*.csv
## ============================================================

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
source(file.path(CODE_V3, "utils_diagnostics.R"))

# Override for 200 draws run
PSI_MAX_DRAWS <<- 200L
RUN_LABEL <<- paste0(RUN_LABEL, "_200draws")
P <- ensure_v3_dirs()


log_time("05_resume", "Starting resume from brms + species trends")

# ── 1. 加载模型与基础数据 ──────────────────────────────────────────────
is_pilot <- Sys.getenv("V3_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL

fit <- safe_read(v3_file("derived", paste0("stMsPGOcc_fit_", run_label), "rds"))
if (is.null(fit)) {
  fit <- safe_read(v3_file("derived", paste0("tMsPGOcc_fit_", run_label), "rds"))
}
if (is.null(fit)) stop("Model fit not found. Run 04 first.")

grid_env_stem <- paste0("grid_environment", GRID_TAG, "_v3")
grid_env <- safe_read(v3_file("derived", grid_env_stem, "rds"))
if (is.null(grid_env)) {
  grid_env <- safe_read(v3_file("derived", "grid_environment_v3", "rds"))
}
if (is.null(grid_env)) {
  grid_env <- safe_read(file.path(DIRS$v2_derived, "grid_environment_dynamic_occupancy.rds"))
}

# ── 2. 提取 psi 后验参数 ──────────────────────────────────────────────
psi_samples <- fit$psi.samples
psi_dim <- dim(psi_samples)
n_draws <- min(dim(psi_samples)[1], PSI_MAX_DRAWS)
draw_idx <- round(seq(1, dim(psi_samples)[1], length.out = n_draws))

n_sp <- dim(psi_samples)[2]
n_sites <- dim(psi_samples)[3]
species <- rownames(fit$y) %||% paste0("sp", seq_len(n_sp))
sites <- colnames(fit$y) %||% paste0("site", seq_len(n_sites))
n_periods <- if (length(psi_dim) >= 4) psi_dim[4] else 1L

message(sprintf("[05_resume] psi_samples dim: [%s]", paste(psi_dim, collapse = "x")))
message(sprintf("[05_resume] n_draws=%d, n_sp=%d, n_sites=%d, n_periods=%d",
                n_draws, n_sp, n_sites, n_periods))

# ── 3. 加载性状与系统发育 ──────────────────────────────────────────────
trait_ext <- safe_read(v3_file("derived", "trait_extended_v3", "rds"))
if (is.null(trait_ext)) {
  warning("[05_resume] trait_extended_v3 not found. Falling back to basic traits.")
  trait_ext <- read_csv_safe(TRAIT_IMPUTED_PATH) |>
    mutate(diet_specialization = NA_real_, habitat_breadth = NA_real_)
}

trait_aligned <- tibble(species = as.character(species)) |>
  left_join(trait_ext %>% mutate(species = as.character(species)), by = "species")

all_trait_vars <- intersect(TRAIT_VARS_ALL, names(trait_aligned))
trait_mat <- prepare_trait_matrix(trait_aligned, all_trait_vars)

phylo <- NULL
if (requireNamespace("ape", quietly = TRUE)) {
  phylo_path <- v3_file("derived", "phylogeny_matched_v3", "rds")
  if (!file.exists(phylo_path)) {
    phylo_path <- file.path(DIRS$v2_derived, "phylogeny_matched.rds")
  }
  phylo <- safe_read(phylo_path)
}

# ── 4. 从 CSV 恢复 trend_summary_df ────────────────────────────────────
trend_csv <- v3_file("results", paste0("table_trend_summary_", run_label), "csv")
if (!file.exists(trend_csv)) {
  stop("trend_summary CSV not found: ", trend_csv,
       "\nPlease ensure the first half of 05 completed successfully.")
}
trend_summary_df <- read_csv(trend_csv, show_col_types = FALSE)
message(sprintf("[05_resume] Loaded trend_summary_df: %d rows", nrow(trend_summary_df)))

# FIX: grid_cell 类型统一（与 05 一致）
if ("grid_cell" %in% names(trend_summary_df)) {
  trend_summary_df <- trend_summary_df |> mutate(grid_cell = as.character(grid_cell))
}
if (!is.null(grid_env) && "grid_cell" %in% names(grid_env)) {
  grid_env$grid_cell <- as.character(grid_env$grid_cell)
}

# ═══════════════════════════════════════════════════════════════════════
# ── 8. brms 驱动回归（加空间 GP）────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════════════
driver_vars <- c("bio4", "bio7", "bio11", "bio13",
                  "elev_mean", "elev_sd", "texture_shannon",
                  "habitat_diversity_shannon",
                  "hfi_mean", "landcover_built", "landcover_cropland",
                  "centroid_lon", "centroid_lat")

for (m in c("corrected_richness", "shannon", "trait_volume", "pd_prob")) {
  message(sprintf("[05_resume] brms driver regression for: %s", m))

  # FIX: 跳过已存在的 brms 结果，避免重复计算
  brms_existing <- v3_file("derived", paste0("brms_driver_trend_", m, "_", run_label), "rds")
  if (file.exists(brms_existing)) {
    message(sprintf("[05_resume] brms result for %s already exists. Skipping.", m))
    next
  }

  trend_m <- trend_summary_df |>
    filter(metric == m) |>
    inner_join(grid_env, by = "grid_cell") |>
    drop_na()

  if (nrow(trend_m) < 30) {
    warning(sprintf("[05_resume] Too few observations for %s (%d), skipping.", m, nrow(trend_m)))
    next
  }

  # 标准化驱动变量
  for (v in driver_vars) {
    if (v %in% names(trend_m)) {
      trend_m[[paste0("z_", v)]] <- scale(trend_m[[v]])[, 1]
    }
  }

  z_vars <- paste0("z_", driver_vars)
  z_vars <- intersect(z_vars, names(trend_m))
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
    warning(sprintf("[05_resume] brms failed for %s: %s", m, e$message))
    NULL
  })

  if (!is.null(brms_fit)) {
    checkpoint_save(brms_fit, paste0("brms_driver_trend_", m, "_", run_label),
                     subdir = "derived")

    # FIX: 用 tryCatch 包裹 check_dharma_gate，防止函数找不到导致崩溃
    dharma_res <- tryCatch(
      check_dharma_gate(
        brms_fit,
        save_dir = DIRS$figures,
        stem = paste0("dharma_driver_trend_", m, "_", run_label)
      ),
      error = function(e) {
        warning(sprintf("[05_resume] DHARMa check failed for %s: %s", m, e$message))
        NULL
      }
    )

    loo_res <- tryCatch(lo(brms_fit), error = function(e) NULL)
    if (!is.null(loo_res)) {
      checkpoint_save(loo_res, paste0("brms_loo_trend_", m, "_", run_label),
                       subdir = "derived")
    }
  }
}

# ── 9. 保存 psi 抽取 ──────────────────────────────────────────────────
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
    "[05_resume] psi_thinned dimension mismatch: got [%s], expected [%s]. Saving anyway.",
    paste(psi_thinned_dim, collapse = "x"),
    paste(expected_dim, collapse = "x")))
}
saveRDS(
  list(psi_samples_thinned = psi_thinned,
       species = species, sites = sites, draw_idx = draw_idx,
       psi_dim = psi_thinned_dim, n_periods = n_periods),
  v3_file("derived", paste0("psi_samples_thinned_", run_label), "rds")
)
message("[05_resume] psi_samples_thinned saved.")

# ═══════════════════════════════════════════════════════════════════════
# ── 10. 物种水平趋势系统分析 ───────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════════════
message("[05_resume] Starting species-level trend analysis (Section 10)")

n_sp_use <- length(species)

sp_trend_arr <- array(NA_real_,
                       dim = c(n_draws, n_sp_use, 2),
                       dimnames = list(draw = NULL, species = species,
                                       method = c("theil_sen", "ols")))

sp_psi_mean <- array(NA_real_,
                      dim = c(n_draws, n_sp_use, n_periods),
                      dimnames = list(draw = NULL, species = species, period = NULL))

message(sprintf("[10a] Computing species trends: %d draws x %d species", n_draws, n_sp_use))

for (d in seq_len(n_draws)) {
  for (sp in seq_len(n_sp_use)) {
    if (length(psi_dim) >= 4) {
      psi_sp <- psi_samples[draw_idx[d], sp, , ]  # sites x periods
    } else {
      psi_sp <- psi_samples[draw_idx[d], sp, ]     # sites（单 period）
      psi_sp <- matrix(psi_sp, ncol = 1)
    }

    psi_mean_t <- colMeans(psi_sp, na.rm = TRUE)
    sp_psi_mean[d, sp, ] <- psi_mean_t

    if (sum(!is.na(psi_mean_t)) >= 3) {
      sp_trend_arr[d, sp, "theil_sen"] <- theil_sen_slope(psi_mean_t)
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

    q <- quantile(vals, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
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

# ── 10c. 物种分组 ──────────────────────────────────────────────────
sp_classify <- sp_trend_df |>
  filter(method == "theil_sen") |>
  mutate(
    trend_class = case_when(
      p_positive > 0.95 ~ "expanding",
      p_negative > 0.95 ~ "contracting",
      TRUE              ~ "stable"
    ),
    trend_grade = case_when(
      p_positive > 0.99 ~ "strong_expanding",
      p_positive > 0.95 ~ "expanding",
      p_negative > 0.99 ~ "strong_contracting",
      p_negative > 0.95 ~ "contracting",
      TRUE              ~ "stable"
    )
  )

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

sp_classify <- sp_classify |>
  left_join(mk_species_df %>% mutate(species = as.character(species)), by = "species")

write_csv(sp_classify, v3_file("results", paste0("table_species_trend_classify_", run_label)))

class_table <- sp_classify |>
  count(trend_class, name = "n_species") |>
  mutate(pct = round(n_species / sum(n_species) * 100, 1))
message(sprintf("[10c] Species classification: %s",
  paste(class_table$trend_class, "=", class_table$n_species, collapse = "; ")))

# ── 10d. 物种趋势 x 性状关联 ───────────────────────────────────────
sp_trait_trend <- sp_classify |>
  select(species, trend_class, trend_grade, mean, p_positive, p_negative, mk_tau, mk_p) |>
  rename(trend_slope = mean) |>
  left_join(trait_ext %>% mutate(species = as.character(species)), by = "species")

write_csv(sp_trait_trend, v3_file("results", paste0("table_species_trend_traits_", run_label)))

trait_corr_results <- list()
for (d in seq_len(min(n_draws, 100))) {
  for (trait_name in all_trait_vars) {
    if (!trait_name %in% names(trait_ext)) next
    trait_vals <- trait_ext[[trait_name]]
    if (all(is.na(trait_vals))) next

    slopes <- sp_trend_arr[d, , "theil_sen"]
    names(slopes) <- species

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

# ── 10e. 物种趋势 x 环境驱动因子 ────────────────────────────────────
if (!is.null(grid_env) && nrow(grid_env) > 0) {
  sp_env_trend <- list()
  psi_final_mean <- apply(psi_samples[draw_idx, , , n_periods], c(2, 3), mean, na.rm = TRUE)

  env_driver_vars <- intersect(
    c("bio4", "bio7", "bio11", "bio13", "elev_mean", "hfi_mean",
      "landcover_built", "landcover_cropland", "centroid_lat", "centroid_lon"),
    names(grid_env)
  )

  for (sp in seq_len(min(n_sp_use, 200))) {
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

  sp_env_trend_df <- sp_classify |>
    select(species, trend_class, trend_slope = mean) |>
    left_join(sp_env_df %>% mutate(species = as.character(species)), by = "species")

  write_csv(sp_env_trend_df, v3_file("results", paste0("table_species_env_trend_", run_label)))

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
expanding_sp  <- sp_classify |> filter(trend_class == "expanding")  |> pull(species)
contracting_sp <- sp_classify |> filter(trend_class == "contracting") |> pull(species)

if (length(expanding_sp) > 0 || length(contracting_sp) > 0) {
  psi_mean_arr <- apply(psi_samples[draw_idx, , , ], c(2, 3, 4), mean, na.rm = TRUE)

  hotspot_df <- list()
  for (t in seq_len(n_periods)) {
    if (length(expanding_sp) > 0) {
      exp_idx <- which(species %in% expanding_sp)
      exp_richness <- colSums(psi_mean_arr[exp_idx, , t], na.rm = TRUE)
    } else {
      exp_richness <- rep(0, n_sites)
    }

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
  message(sprintf("[10f] Hotspot map data: %d grids x %d periods", n_sites, n_periods))
}

message(sprintf("[10] Species-level analysis complete: %d expanding, %d stable, %d contracting",
  sum(sp_classify$trend_class == "expanding"),
  sum(sp_classify$trend_class == "stable"),
  sum(sp_classify$trend_class == "contracting")))

log_time("05_resume", "DONE")
