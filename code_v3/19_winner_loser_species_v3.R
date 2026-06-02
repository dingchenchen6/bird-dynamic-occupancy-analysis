#!/usr/bin/env Rscript
## 19_winner_loser_species_v3.R  —  v3 物种赢家/输家
## 基于已生成的 brms 物种系数 / spOccupancy beta.samples，识别 winner/loser
## 输出：figures_v3/fig_v3_winner_loser_* + results_v3/table_species_winners_losers_v3.csv

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble); library(stringr)
  library(forcats); library(ggplot2); library(patchwork)
  library(ggridges); library(ggdist); library(ggbeeswarm)
})

CODE_V3 <- Sys.getenv("V3_CODE_DIR", "/home/dingchenchen/bird_dynamic_occupancy_analysis/code_v3")
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_mapping.R"))
ensure_v3_dirs()
log_time("19", "winner/loser species analysis")

# ---- 1. species year_scaled 系数 ----
RUN_LABEL <- Sys.getenv("V3_RUN_LABEL", RUN_LABEL)
coef_csv <- v3_file("results", paste0("table_species_coefficients_", RUN_LABEL), "csv")
if (!file.exists(coef_csv)) {
  v2 <- file.path(DIRS$v2_results, "table_species_coefficients_v2_full_200sp_ar1.csv")
  if (file.exists(v2)) coef_csv <- v2
  else stop("Cannot find species coefficient table")
}
species_coef <- read_csv(coef_csv, show_col_types = FALSE)
trait_ext <- safe_read(v3_file("derived", "trait_extended_v3", "rds"))
if (is.null(trait_ext)) trait_ext <- safe_read(file.path(DIRS$v2_derived, "trait_extended.rds"))

# 兼容不同列名：mean/median, l95/q025, u95/q975
norm_cols <- function(df){
  rename_map <- c(mean="mean", median="median", l95="l95", u95="u95",
                  q025="l95", q975="u95", sig95="sig95")
  for (old in names(rename_map)) if (old %in% names(df)) names(df)[names(df)==old] <- rename_map[[old]]
  df
}
species_coef <- norm_cols(species_coef)
yr <- species_coef |> filter(covariate == "year_scaled") |>
  transmute(species,
            year_slope = if ("mean" %in% names(species_coef)) mean else median,
            year_l95 = l95, year_u95 = u95,
            year_sig = if ("sig95" %in% names(species_coef)) sig95 else (l95>0 | u95<0))

ws <- yr |> left_join(trait_ext, by = "species") |>
  mutate(status = factor(case_when(
    year_sig & year_slope > 0 ~ "Winner",
    year_sig & year_slope < 0 ~ "Loser",
    TRUE ~ "Stable"), levels = c("Winner","Stable","Loser")))
log_time("19", sprintf("Winners=%d Stable=%d Losers=%d",
  sum(ws$status=="Winner"), sum(ws$status=="Stable"), sum(ws$status=="Loser")))
write_csv(ws, v3_file("results", "table_species_winners_losers_v3", "csv"))

pal <- c(Winner="#8B2E1E", Loser="#0E5A78", Stable="grey55")

# ---- 2. 森林图 top60 ----
top_w <- ws |> filter(status=="Winner") |> arrange(desc(year_slope)) |> head(30)
top_l <- ws |> filter(status=="Loser") |> arrange(year_slope) |> head(30)
top <- bind_rows(top_w, top_l) |> mutate(species = forcats::fct_reorder(species, year_slope))
p_forest <- ggplot(top, aes(year_slope, species, colour=status)) +
  geom_vline(xintercept=0, linetype=2, colour="grey55") +
  geom_errorbarh(aes(xmin=year_l95, xmax=year_u95), height=0, linewidth=0.5, alpha=0.85) +
  geom_point(size=2.2) +
  scale_colour_manual(values=pal, guide="none") +
  labs(title="Top 30 winners and top 30 losers",
       subtitle=sprintf("Per-species year_scaled slope | %s", RUN_LABEL),
       x="Year-trend coefficient (logit)", y=NULL) +
  theme_nature_pub(7.5) + theme(axis.text.y=element_text(size=6))
save_nature(p_forest, "fig_v3_winner_loser_forest", width_mm=140, height_mm=200)

# ---- 3. 雨林 by trait group ----
group_cols <- intersect(c("Trophic.Niche","Habitat","Primary.Lifestyle"), names(ws))
plots_g <- lapply(group_cols, function(g){
  d <- ws |> filter(!is.na(.data[[g]])) |>
    mutate(grp = forcats::fct_reorder(.data[[g]], year_slope, .fun=median))
  ggplot(d, aes(year_slope, grp, fill=status, colour=status)) +
    geom_vline(xintercept=0, linetype=2, colour="grey55") +
    ggdist::stat_halfeye(thickness=0.55, .width=c(0.5,0.95), slab_alpha=0.45,
                          slab_colour=NA, side="right", justification=-0.15,
                          show.legend=FALSE) +
    ggbeeswarm::geom_quasirandom(groupOnX=FALSE, alpha=0.45, size=0.6, width=0.18,
                                  show.legend=FALSE) +
    scale_colour_manual(values=pal) + scale_fill_manual(values=pal) +
    labs(title=g, x="Year-trend coefficient", y=NULL) + theme_nature_pub(7) +
    theme(legend.position="top")
})
if (length(plots_g) > 0) {
  p_rain <- patchwork::wrap_plots(plots_g, ncol=1) +
    patchwork::plot_annotation(title="Winner/loser by trait group", theme=theme_nature_pub(8))
  save_nature(p_rain, "fig_v3_winner_loser_raincloud_by_trait", width_mm=183, height_mm=240)
}

# ---- 4. 山脊 by family ----
fam_col <- intersect(c("family_lat","family"), names(ws))[1]
if (!is.na(fam_col)) {
  fam_top <- ws |> filter(!is.na(.data[[fam_col]])) |>
    count(.data[[fam_col]], sort=TRUE) |> filter(n>=4) |> head(15) |>
    pull(1)
  d_fam <- ws |> filter(.data[[fam_col]] %in% fam_top) |>
    mutate(fam = forcats::fct_reorder(.data[[fam_col]], year_slope, .fun=median))
  p_ridge <- ggplot(d_fam, aes(year_slope, fam, fill=after_stat(x))) +
    geom_vline(xintercept=0, linetype=2, colour="grey55") +
    ggridges::geom_density_ridges_gradient(scale=2.0, rel_min_height=0.008,
                                            colour="white", linewidth=0.2, alpha=0.95) +
    scale_fill_gradient2(low="#0E5A78", mid="#F2E8D8", high="#8B2E1E", midpoint=0, guide="none") +
    labs(title="Year-trend distributions by family",
         subtitle=sprintf("Top 15 families (>=4 species each) | %s", RUN_LABEL),
         x="Year-trend coefficient", y=NULL) + theme_nature_pub(7.5)
  save_nature(p_ridge, "fig_v3_winner_loser_ridge_by_family", width_mm=140, height_mm=160)
}

# ---- 5. trait scatter ----
tr_cols <- intersect(c("body_mass_g","clutch_size","longevity_y","avonet_hwi",
                        "avonet_range_size","Habitat.Density","habitat_breadth",
                        "diet_specialization"), names(ws))
if (length(tr_cols) > 0) {
  sl <- ws |> select(year_slope, status, all_of(tr_cols)) |>
    pivot_longer(all_of(tr_cols), names_to="trait", values_to="val") |>
    filter(is.finite(val))
  p_sc <- ggplot(sl, aes(val, year_slope)) +
    geom_hline(yintercept=0, linetype=2, colour="grey55") +
    geom_point(aes(colour=status), alpha=0.5, size=0.9) +
    geom_smooth(method="lm", se=TRUE, colour="#2A2A2A",
                fill="#0E5A78", alpha=0.12, linewidth=0.5) +
    facet_wrap(~trait, scales="free_x", ncol=3) +
    scale_colour_manual(values=pal, name=NULL) +
    labs(title="Species traits vs occupancy time trend",
         subtitle=sprintf("Run=%s", RUN_LABEL),
         x=NULL, y="Year-trend coefficient") +
    theme_nature_pub(7.5) + theme(legend.position="top")
  save_nature(p_sc, "fig_v3_trait_vs_trend_scatter", width_mm=183, height_mm=160)
}

log_time("19", "done")
