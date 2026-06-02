#!/usr/bin/env Rscript
## 14_species_trait_regression.R  —  v3 Q4 性状回归
##
## 物种水平占有率趋势的性状解释力
## v3 新增：diet_specialization + habitat_breadth
## brms 生态位模型含系统发育随机效应

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(brms); library(cmdstanr)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
P <- ensure_v3_dirs()

log_time("14", "Starting species-trait regression (Q4)")

# ── 1. 加载种级趋势 ──────────────────────────────────────────────────
is_pilot <- Sys.getenv("V3_PILOT", "0") == "1"
run_label <- if (is_pilot) PILOT_LABEL else RUN_LABEL

trend_df <- read_csv_safe(v3_file("results",
                                   paste0("table_species_trend_", run_label, "_extended"))) |>

  filter(method == "theil_sen") |> select(species, trend_i = mean, trend_q025 = q025, trend_q975 = q975)

# ── 2. 加载扩展性状 ──────────────────────────────────────────────────
trait_ext <- safe_read(v3_file("derived", "trait_extended_v3", "rds"))
if (is.null(trait_ext)) {
  stop("trait_extended_v3 not found. Run 03b first.")
}

# 合并
reg_df <- trend_df |>
  left_join(trait_ext, by = "species") |>
  drop_na(trend_i)

# 标准化性状
trait_cols_reg <- c("body_mass_g", "clutch_size", "avonet_hwi",
                     "avonet_range_size",
                     "diet_specialization", "habitat_breadth")
trait_cols_reg <- intersect(trait_cols_reg, names(reg_df))

# FIX #14: 使用 log10 变换（与 utils_diversity.R 的 prepare_trait_matrix 一致）
# 旧版 log1p 与 diversity 计算不兼容
log_vars_14 <- intersect(TRAIT_VARS_LOG10, trait_cols_reg)
no_log_vars_14 <- intersect(TRAIT_VARS_NO_LOG, trait_cols_reg)
# 其余未在 TRAIT_VARS_LOG10 / TRAIT_VARS_NO_LOG 中的变量不变换
other_vars_14 <- setdiff(trait_cols_reg, c(log_vars_14, no_log_vars_14))

for (cc in log_vars_14) {
  reg_df[[paste0("z_", cc)]] <- scale(if_else(reg_df[[cc]] > 0, log10(reg_df[[cc]]), NA_real_))[, 1]
}
for (cc in no_log_vars_14) {
  reg_df[[paste0("z_", cc)]] <- scale(reg_df[[cc]])[, 1]
}

z_trait_cols <- paste0("z_", trait_cols_reg)
z_trait_cols <- intersect(z_trait_cols, names(reg_df))

message(sprintf("[14] %d species with complete trait data (%d traits)",
                nrow(reg_df), length(z_trait_cols)))

# ── 3. 加载系统发育 ──────────────────────────────────────────────────
# v3 命名一致性 fallback:phylogeny_matched_v3 → phylogeny_mctavish_matched_v3 → v2 旧路径
phylo <- safe_read(v3_file("derived", "phylogeny_matched_v3", "rds"))
if (is.null(phylo)) {
  phylo <- safe_read(v3_file("derived", "phylogeny_mctavish_matched_v3", "rds"))
}
if (is.null(phylo)) {
  phylo <- safe_read(file.path(DIRS$v2_derived, "phylogeny_matched.rds"))
}

# 系统发育距离矩阵
A <- NULL
if (!is.null(phylo)) {
  matched_sp <- intersect(reg_df$species, phylo$tip.label)
  if (length(matched_sp) >= 30) {
    pruned <- ape::drop.tip(phylo, setdiff(phylo$tip.label, matched_sp))
    A <- ape::vcv.phylo(pruned, corr = TRUE)
    reg_df <- reg_df |> filter(species %in% matched_sp)
    message(sprintf("[14] Phylogeny: %d/%d species matched", length(matched_sp), nrow(reg_df)))
  }
}

# ── 4. brms 性状回归 ──────────────────────────────────────────────────
# 公式：trend_i ~ z_body_mass + z_hwi + z_range_size + z_clutch_size
#              + z_diet_specialization + z_habitat_breadth
#              + (1 | gr(species, cov = A))

# Build dynamic formula from available traits
trait_terms <- paste(z_trait_cols, collapse = " + ")
if (trait_terms == "") stop("No trait columns available for regression")

if (!is.null(A)) {
  brms_formula <- bf(as.formula(paste("trend_i ~", trait_terms, "+ (1 | gr(species, cov = A))")))
} else {
  brms_formula <- bf(as.formula(paste("trend_i ~", trait_terms)))
  warning("[14] Phylogeny not available. Fitting without phylogenetic random effect.")
}

brms_prior <- if (!is.null(A)) {
  c(
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(normal(0, 1), class = "Intercept")
  )
} else {
  c(
    prior(normal(0, 1), class = "b"),
    prior(normal(0, 1), class = "Intercept")
  )
}

fit_14 <- tryCatch({
  brm(
    formula    = brms_formula,
    data       = reg_df,
    data2      = if (!is.null(A)) list(A = A) else NULL,
    family     = gaussian(),
    prior      = brms_prior,
    iter       = BRMS_ITER,
    warmup     = BRMS_WARMUP,
    chains     = BRMS_CHAINS,
    cores      = BRMS_CHAINS,
    control    = list(adapt_delta = BRMS_ADAPT_DELTA,
                       max_treedepth = BRMS_MAX_TREED),
    seed       = BRMS_SEED,
    backend    = "cmdstanr",
    silent     = 2
  )
}, error = function(e) {
  warning(sprintf("[14] brms trait regression failed: %s", e$message))
  NULL
})

if (!is.null(fit_14)) {
  checkpoint_save(fit_14, paste0("brms_trait_regression_", run_label), subdir = "derived")

  # 系数摘要
  coef_df <- posterior_summary(fit_14, pars = "^b_") |>
    as.data.frame() |>
    tibble::rownames_to_column("term") |>
    as_tibble() |>
    rename(estimate = Estimate, q025 = `Q2.5`, q975 = `Q97.5`) |>
    mutate(term = gsub("^b_", "", term))

  write_csv(coef_df, v3_file("results", paste0("table_trait_regression_coefs_", run_label)))

  # R2
  r2 <- tryCatch(bayes_R2(fit_14), error = function(e) NULL)
  if (!is.null(r2) && is.matrix(r2)) {
    cn <- colnames(r2)
    if ("R2_Conditional" %in% cn) {
      message(sprintf("[14] R2_conditional = %.3f, R2_marginal = %.3f",
                       mean(r2[, "R2_Conditional"], na.rm = TRUE),
                       mean(r2[, "R2_Marginal"], na.rm = TRUE)))
    } else if ("Estimate" %in% cn) {
      message(sprintf("[14] R2 = %.3f [%.3f, %.3f]",
                       r2["R2", "Estimate"], r2["R2", "Q2.5"], r2["R2", "Q97.5"]))
    } else {
      message("[14] R2 matrix format unexpected")
    }
  } else if (!is.null(r2) && is.numeric(r2)) {
    message(sprintf("[14] R2 = %.3f", mean(r2, na.rm = TRUE)))
  } else {
    message("[14] Could not compute R2")
  }

  # LOO
  loo_res <- tryCatch(brms::loo(fit_14), error = function(e) NULL)
  if (!is.null(loo_res)) {
    message(sprintf("[14] LOOIC = %.1f", loo_res$estimates["looic", "Estimate"]))
  }

  # DHARMa
# check_dharma_gate(fit_14, save_dir = DIRS$figures,
  #                    stem = paste0("dharma_trait_regression_", run_label))
}

log_time("14", "DONE")
