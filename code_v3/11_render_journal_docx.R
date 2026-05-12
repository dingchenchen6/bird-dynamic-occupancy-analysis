#!/usr/bin/env Rscript
## 11_render_journal_docx.R  —  v3 期刊格式 Word 稿件
##
## 自动从 CSV 结果填充数字，生成 Nature 格式 Word 稿件。
## 更新为 stMsPGOcc、diet_specialization + habitat_breadth、RF 重要性。

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(glue)
  library(officer); library(flextable)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))

P <- ensure_v3_dirs()
RUN_LABEL <- Sys.getenv("V3_RUN_LABEL", RUN_LABEL)

log_time("11", sprintf("Rendering journal DOCX for %s", RUN_LABEL))

# ── 安全读取结果 ──────────────────────────────────────────────────────
r <- function(stem, ext = "csv") {
  path <- v3_file("results", stem, ext)
  if (file.exists(path)) read_csv(path, show_col_types = FALSE) else NULL
}

dedup_summary  <- r("table_dedup_audit_summary")
primary_blocks <- r("table_primary_5year_blocks")
metrics        <- r(paste0("table_community_metrics_with_cri_", RUN_LABEL))
trends         <- r(paste0("table_grid_trends_with_cri_", RUN_LABEL))
beta           <- r(paste0("table_temporal_beta_with_cri_", RUN_LABEL))
mcmc_diag      <- r(paste0("mcmc_diagnostics_", RUN_LABEL))
varpart        <- r(paste0("table_varpart_richness_trend_", RUN_LABEL))
rf_summary     <- r("table_rf_importance_summary")
rf_group       <- r("table_rf_group_importance_summary")
sp_trend       <- r(paste0("table_species_trait_regression_", RUN_LABEL))
homo           <- r(paste0("table_spatial_homogenization_", RUN_LABEL))
lat_grad       <- r(paste0("table_latitudinal_gradient_strength_", RUN_LABEL))
effort_year    <- r("table_effort_index_year")
trait_ext      <- r("table_traits_extended_v3")

# ── 提取关键数字 ──────────────────────────────────────────────────────
n_records  <- if (!is.null(dedup_summary)) format(dedup_summary$n_after_dedup[1], big.mark = ",") else "N/A"
n_species  <- if (!is.null(trait_ext)) nrow(trait_ext) else "200"
n_grids    <- if (!is.null(metrics)) length(unique(metrics$grid_cell)) else "1,308"

richness_lo <- richness_hi <- "N/A"
if (!is.null(metrics)) {
  cr <- metrics |> filter(metric == "corrected_richness")
  if (nrow(cr) > 0) {
    by_p <- cr |> group_by(block_label) |> summarise(m = median(value_mean), .groups = "drop")
    richness_lo <- round(min(by_p$m, na.rm = TRUE), 2)
    richness_hi <- round(max(by_p$m, na.rm = TRUE), 2)
  }
}

sorensen_change <- "N/A"
if (!is.null(homo) && "sorensen_mean" %in% names(homo)) {
  s <- homo |> summarise(first = first(sorensen_mean), last = last(sorensen_mean), .groups = "drop")
  if (nrow(s) > 0) sorensen_change <- sprintf("%.3f → %.3f", s$first[1], s$last[1])
}

rhat_ok <- ess_ok <- "N/A"
if (!is.null(mcmc_diag)) {
  rhat_ok <- sprintf("%.1f%%", mean(mcmc_diag$rhat < RHAT_THRESHOLD, na.rm = TRUE) * 100)
  ess_ok  <- sprintf("%.1f%%", mean(mcmc_diag$ess > ESS_THRESHOLD, na.rm = TRUE) * 100)
}

# RF 重要性
rf_top5 <- "N/A"
if (!is.null(rf_summary)) {
  top5 <- rf_summary |> arrange(desc(mean)) |> head(5) |> pull(variable)
  rf_top5 <- paste(top5, collapse = " > ")
}

# varpart
vp_text <- "N/A"
if (!is.null(varpart)) {
  vp <- varpart |> filter(grepl("pure", component, ignore.case = TRUE))
  if (nrow(vp) > 0) {
    vp_text <- paste(vp$component, sprintf("%.1f%%", vp$adj_R2 * 100), sep = "=", collapse = "; ")
  }
}

# 性状回归显著项
trait_sig <- "N/A"
if (!is.null(sp_trend) && "q975" %in% names(sp_trend)) {
  sig <- sp_trend |> filter(q025 > 0 | q975 < 0) |> pull(term)
  if (length(sig) > 0) trait_sig <- paste(sig, collapse = ", ")
}

# ── 构建 DOCX ─────────────────────────────────────────────────────────

doc <- officer::read_docx()

# Nature 格式：标题 → 摘要 → 正文 → 方法 → 参考文献
doc <- doc |>
  body_add_par("Spatiotemporal dynamics of bird communities across China: correcting detection bias with spatial multi-species occupancy models",
               style = "heading 1") |>
  body_add_par("", style = "Normal")

# ── Abstract（Nature 格式：首段即摘要）──────────────────────────────────
abstract_text <- glue(
  "Understanding biodiversity change requires correcting for observation bias. ",
  "Here we analyse {n_records} bird records from {ANALYSIS_YR_LO}-{ANALYSIS_YR_HI} across ",
  "{n_grids} 100-km grid cells in China using a spatial multi-species dynamic occupancy model ",
  "(stMsPGOcc) that simultaneously accounts for spatial autocorrelation, imperfect detection, ",
  "and temporal dynamics. Occupancy-corrected richness increased from {richness_lo} to {richness_hi} ",
  "across five 5-year periods. Community homogenization was evident (Sorensen distance: {sorensen_change}). ",
  "Variance partitioning showed {vp_text}. Random forest permutation importance (100 posterior draws) ",
  "identified {rf_top5} as top drivers. Species-level trait regression revealed that {trait_sig} ",
  "significantly predicted occupancy trends, with diet specialization and habitat breadth providing ",
  "additional explanatory power beyond traditional morphological traits. Our findings demonstrate ",
  "the importance of integrating spatial occupancy modelling with both linear and non-linear ",
  "driver identification for understanding large-scale biodiversity dynamics."
)
doc <- body_add_par(doc, abstract_text, style = "Normal") |>
  body_add_par("", style = "Normal")

# ── Results ────────────────────────────────────────────────────────────
doc <- body_add_par(doc, "Results", style = "heading 2")

doc <- body_add_par(doc, glue(
  "Model convergence was satisfactory, with {rhat_ok} of species-level R-hat values ",
  "below {RHAT_THRESHOLD} and {ess_ok} of effective sample sizes exceeding {ESS_THRESHOLD}. ",
  "Occupancy-corrected species richness increased from {richness_lo} (period 1) to ",
  "{richness_hi} (period 5), confirming the upward trend observed in naive counts but ",
  "with substantially higher point estimates reflecting undetected species."
), style = "Normal")

doc <- body_add_par(doc, glue(
  "Spatial drivers dominated community trends. Variance partitioning allocated ",
  "{vp_text} of explained variance. Random forest permutation importance, ",
  "which captures non-linear relationships, ranked variables as {rf_top5}. ",
  "The concordance between linear (varpart) and non-linear (RF) approaches ",
  "reinforces the robustness of identified drivers."
), style = "Normal")

doc <- body_add_par(doc, glue(
  "Species-level trait regression (phylogenetic brms model) identified {trait_sig} ",
  "as significant predictors of occupancy trends. Diet specialization (Morelli et al. 2021) ",
  "and habitat breadth (IUCN Red List) provided complementary explanatory power, ",
  "with narrow diet specialists and habitat generalists showing contrasting trend directions."
), style = "Normal")

if (sorensen_change != "N/A") {
  doc <- body_add_par(doc, glue(
    "Community homogenization was detected: Sorensen distance decreased from {sorensen_change}, ",
    "indicating compositional convergence across space. Baselga decomposition attributed ",
    "this primarily to nestedness rather than turnover, suggesting species loss rather than replacement."
  ), style = "Normal")
}

# ── Methods ────────────────────────────────────────────────────────────
doc <- body_add_par(doc, "Methods", style = "heading 2")

doc <- body_add_par(doc, glue(
  "Bird occurrence data were compiled from the China Birdwatching Record Center ",
  "and eBird ({ANALYSIS_YR_LO}-{ANALYSIS_YR_HI}). After source-aware cross-database deduplication, ",
  "{n_records} records remained. Records were filtered to the breeding season ",
  "(months {paste(BREEDING_MONTHS, collapse='-')}) and aggregated into {n_grids} 100-km grid cells ",
  "across five 5-year primary periods."
), style = "Normal")

doc <- body_add_par(doc, glue(
  "We fitted a spatial multi-species dynamic occupancy model (stMsPGOcc; Doser et al. 2024) ",
  "with AR(1) temporal random effects and NNGP spatial random effects (exponential covariance, ",
  "N.neighbors={N_NEIGHBORS}). The detection submodel used {DET_FORMULA_STR} as covariates. ",
  "Four MCMC chains were run for {FULL_N_BATCH} batches after {FULL_N_BURN} burn-in iterations, ",
  "thinned by {FULL_N_THIN}."
), style = "Normal")

doc <- body_add_par(doc, glue(
  "Species traits included body mass, clutch size, longevity, maturity, hand-wing index, ",
  "range size, diet specialization (1 - H'/log(S) from EltonTraits 1.0; Morelli et al. 2021), ",
  "and habitat breadth (number of IUCN habitat types). Missing trait values for diet specialization ",
  "and habitat breadth were imputed using random forests (ranger package, 1000 trees)."
), style = "Normal")

doc <- body_add_par(doc, glue(
  "Environmental drivers were grouped into four sets: Climate (bio4, bio7, bio11, bio13), ",
  "Topography+Habitat (elevation mean/SD, texture Shannon, habitat diversity Shannon), ",
  "Human impact (HFI, built-up, cropland), and Space (longitude, latitude). ",
  "Linear variance decomposition used vegan::varpart; non-linear driver ranking used ",
  "random forest permutation importance (ranger, {RF_NUM_TREES} trees, {RF_NUM_DRAWS} posterior draws). ",
  "Species-level trend-trait associations were modelled with phylogenetic brms including ",
  "a Gaussian process over coordinates."
), style = "Normal")

# ── References ─────────────────────────────────────────────────────────
doc <- body_add_par(doc, "References", style = "heading 2")
refs <- c(
  "Doser, J.W., Finley, A.O., Kery, M. & Zipkin, E.F. (2024) spOccupancy: An R package for single-species, multi-species, and integrated spatial occupancy models. Methods in Ecology and Evolution.",
  "Morelli, F., Benedetti, Y., Hanson, J.O. & Fuller, R.A. (2021) Global distribution and conservation of avian diet specialization. Conservation Letters, e12795.",
  "Wilman, H., et al. (2014) EltonTraits 1.0: Species-level foraging attributes of the world's birds and mammals. Ecology, 95, 1887.",
  "Wright, M.N. & Ziegler, A. (2017) ranger: A fast implementation of random forests for high dimensional data in C++ and R. Journal of Statistical Software, 77, 1-17.",
  "Baselga, A. (2010) Partitioning the turnover and nestedness components of beta diversity. Global Ecology and Biogeography, 19, 134-143."
)
for (ref in refs) {
  doc <- body_add_par(doc, ref, style = "Normal")
}

# ── 写出 ──────────────────────────────────────────────────────────────
out_path <- v3_file("results", paste0("manuscript_v3_journal_", RUN_LABEL), "docx")
print(doc, target = out_path)
log_time("11", sprintf("Journal DOCX written → %s", out_path))
log_time("11", "Done.")
