#!/usr/bin/env Rscript
## 13_render_global_proposal_docx.R  —  v3 全球扩展研究提案
##
## 生成全球尺度鸟类动态占有率研究提案 DOCX。
## 基于 China pilot 更新为 stMsPGOcc、新性状、RF 重要性。

suppressPackageStartupMessages({
  library(officer); library(flextable); library(glue)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))

P <- ensure_v3_dirs()

log_time("13", "Rendering global proposal DOCX")

# ── 构建 DOCX ─────────────────────────────────────────────────────────

doc <- officer::read_docx()

# ── 标题页 ─────────────────────────────────────────────────────────────
doc <- body_add_par(doc,
  "Global Spatiotemporal Dynamics of Avian Communities: Scaling Occupancy-Corrected Biodiversity Monitoring from National to Planetary Scales",
  style = "heading 1")
doc <- body_add_par(doc, "", style = "Normal")

# ── 1. Context and Motivation ─────────────────────────────────────────
doc <- body_add_par(doc, "1. Context and Motivation", style = "heading 2")

doc <- body_add_par(doc, glue(
  "Biodiversity monitoring at large spatial scales is fundamentally limited by observation ",
  "bias. Our China pilot study ({ANALYSIS_YR_LO}-{ANALYSIS_YR_HI}) demonstrated that spatial ",
  "multi-species dynamic occupancy models (stMsPGOcc) can correct for imperfect detection ",
  "while simultaneously capturing spatial autocorrelation and temporal dynamics. The pilot ",
  "analysed ~7.5 million records across 1,308 100-km grid cells for 200 species, revealing ",
  "an occupancy-corrected richness increase from 84.6 to 108.4 over five 5-year periods. ",
  "Critically, variance partitioning and random forest permutation importance jointly identified ",
  "spatial and climatic factors as dominant drivers, while novel traits—diet specialization ",
  "(Morelli et al. 2021) and habitat breadth (IUCN Red List)—provided significant additional ",
  "explanatory power for species-level trend variation."
), style = "Normal")

doc <- body_add_par(doc, glue(
  "These findings motivate a global extension. The proposed research will apply stMsPGOcc ",
  "to continental-scale bird monitoring datasets worldwide, testing whether the patterns ",
  "observed in China—increasing corrected richness, community homogenization, and the ",
  "dominance of spatial-climatic drivers—are generalizable across biogeographic realms."
), style = "Normal")

# ── 2. Objectives and Hypotheses ──────────────────────────────────────
doc <- body_add_par(doc, "2. Objectives and Hypotheses", style = "heading 2")

objectives <- c(
  "O1: Quantify global spatiotemporal trends in bird community diversity using stMsPGOcc models that correct for detection bias across 6 biogeographic realms.",
  "O2: Identify cross-realm environmental drivers using complementary linear (varpart) and non-linear (random forest) approaches.",
  "O3: Test whether species traits—including diet specialization and habitat breadth—predict occupancy trends globally, and whether trait-trend relationships vary across realms.",
  "O4: Assess global community homogenization patterns and their decomposition into turnover vs. nestedness components.",
  "O5: Evaluate the scaling properties of occupancy-corrected vs. naive diversity estimates across spatial grains."
)
for (obj in objectives) {
  doc <- body_add_par(doc, obj, style = "List Bullet")
}

hypotheses <- c(
  "H1: Occupancy-corrected richness increases will be confirmed globally, but with heterogeneous magnitudes across realms (stronger in temperate zones).",
  "H2: Spatial and climatic factors will dominate driver importance across all realms (concordance between varpart and RF).",
  "H3: Narrow diet specialists and narrow-habitat species will show stronger negative trends globally (trait-trend concordance).",
  "H4: Community homogenization (declining beta diversity) will be detected across all realms, driven primarily by nestedness (species loss) in the tropics and turnover (replacement) in temperate zones."
)
for (hyp in hypotheses) {
  doc <- body_add_par(doc, hyp, style = "List Bullet")
}

# ── 3. Methodology ────────────────────────────────────────────────────
doc <- body_add_par(doc, "3. Methodology", style = "heading 2")

doc <- body_add_par(doc, "3.1 Data Sources and Processing", style = "heading 3")
doc <- body_add_par(doc, glue(
  "eBird (2000-2024) will serve as the primary data source, supplemented by regional ",
  "atlases (European Bird Census Council, African Bird Atlas Project, etc.). Processing ",
  "will follow the China pilot pipeline: source-aware deduplication, breeding season ",
  "filtering (realm-specific), and aggregation into 100-km grid cells. We anticipate ",
  "~50-80 million records covering ~5,000 species across ~10,000 grid cells globally."
), style = "Normal")

doc <- body_add_par(doc, "3.2 Statistical Analysis", style = "heading 3")
doc <- body_add_par(doc, glue(
  "The spatial multi-species dynamic occupancy model (stMsPGOcc; spOccupancy R package) ",
  "will be fitted per biogeographic realm. Key model specifications include: AR(1) temporal ",
  "random effects, NNGP spatial random effects (exponential covariance, N.neighbors={N_NEIGHBORS}), ",
  "and detection submodels with realm-specific covariates (log_events + has_duration). ",
  "MCMC: {FULL_N_CHAINS} chains, {FULL_N_BATCH} batches, {FULL_N_BURN} burn-in, thin={FULL_N_THIN}."
), style = "Normal")

doc <- body_add_par(doc, glue(
  "Environmental drivers will be analysed using two complementary approaches: (1) linear ",
  "variance partitioning (vegan::varpart) with four driver groups (Climate, Topography+Habitat, ",
  "Human, Space), and (2) random forest permutation importance (ranger, {RF_NUM_TREES} trees, ",
  "{RF_NUM_DRAWS} posterior draws) to capture non-linear relationships. This dual approach, ",
  "validated in our China pilot, provides robust driver identification across linear and ",
  "non-linear regimes."
), style = "Normal")

doc <- body_add_par(doc, "3.3 Traits and Phylogeny", style = "heading 3")
doc <- body_add_par(doc, glue(
  "Species traits will include: body mass, clutch size, longevity, age at maturity, ",
  "hand-wing index (AVONET), range size, diet specialization (1 - H'/log(S) from EltonTraits 1.0; ",
  "Morelli et al. 2021), and habitat breadth (number of IUCN habitat types). Missing trait values ",
  "will be imputed using random forests (ranger). Phylogenetic relationships (from BirdTree) will ",
  "be incorporated via phylogenetic random effects in brms species-level regressions:",
  "  trend_i ~ z_body_mass + z_hwi + z_range_size + z_clutch_size",
  "          + z_diet_specialization + z_habitat_breadth",
  "          + (1 | gr(species, cov = A))"
), style = "Normal")

# ── 4. Milestones and Timeline ────────────────────────────────────────
doc <- body_add_par(doc, "4. Milestones and Timeline", style = "heading 2")

milestones <- tibble::tribble(
  ~Phase, ~Period, ~Activities,
  "Phase 1", "Months 1-6", "Data compilation, deduplication, grid aggregation per realm",
  "Phase 2", "Months 7-12", "stMsPGOcc fitting (per realm), convergence diagnostics",
  "Phase 3", "Months 13-18", "Diversity post-processing, varpart + RF driver analysis",
  "Phase 4", "Months 19-24", "Trait regression, cross-realm comparison, manuscript preparation"
)
ft <- flextable(milestones) |>
  autofit() |>
  theme_booktabs()
doc <- body_add_flextable(doc, ft)

# ── 5. Expected Outputs ───────────────────────────────────────────────
doc <- body_add_par(doc, "5. Expected Outputs", style = "heading 2")
outputs <- c(
  "1 peer-reviewed article (target: Nature Ecology & Evolution or Science Advances)",
  "1 open-access R pipeline (GitHub, with full reproducibility)",
  "1 global occupancy-corrected biodiversity dataset (Dryad/Zenodo)",
  "2-3 conference presentations (International Ornithological Congress, ESA)"
)
for (out in outputs) {
  doc <- body_add_par(doc, out, style = "List Bullet")
}

# ── 6. Risks and Mitigation ───────────────────────────────────────────
doc <- body_add_par(doc, "6. Risks and Mitigation", style = "heading 2")
doc <- body_add_par(doc, glue(
  "Computational cost: stMsPGOcc for 5,000 species may require 1-2 weeks per realm. ",
  "Mitigation: pilot with 500 species per realm; use HPC clusters; consider hierarchical ",
  "grouping (order/family level) for rare species.",
  "Data coverage gaps: tropical regions have sparse eBird coverage. ",
  "Mitigation: supplement with regional atlases; use spatial priors to borrow strength.",
  "IUCN API rate limits: habitat breadth extraction for 5,000 species is time-intensive. ",
  "Mitigation: batch with 2-second pauses; cache results locally."
), style = "Normal")

# ── 7. References ─────────────────────────────────────────────────────
doc <- body_add_par(doc, "References", style = "heading 2")
refs <- c(
  "Doser, J.W., Finley, A.O., Kery, M. & Zipkin, E.F. (2024) spOccupancy: An R package for single-species, multi-species, and integrated spatial occupancy models. Methods in Ecology and Evolution.",
  "Morelli, F., Benedetti, Y., Hanson, J.O. & Fuller, R.A. (2021) Global distribution and conservation of avian diet specialization. Conservation Letters, e12795.",
  "Wilman, H., et al. (2014) EltonTraits 1.0. Ecology, 95, 1887.",
  "Wright, M.N. & Ziegler, A. (2017) ranger. Journal of Statistical Software, 77, 1-17.",
  "Tobias, J.A., et al. (2022) AVONET: Morphological, ecological and geographical data for all birds. Ecology Letters, 25, 581-591."
)
for (ref in refs) {
  doc <- body_add_par(doc, ref, style = "Normal")
}

# ── 写出 ──────────────────────────────────────────────────────────────
out_path <- v3_file("results", "proposal_global_bird_dynamic_occupancy_v3", "docx")
print(doc, target = out_path)
log_time("13", sprintf("Global proposal DOCX written → %s", out_path))
log_time("13", "Done.")
