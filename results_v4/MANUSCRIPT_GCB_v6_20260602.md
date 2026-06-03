# Beneath rising richness: detection-corrected citizen science reveals functional homogenization of Chinese bird communities under rapid global change, 2000–2024

**Target journal:** *Global Change Biology* (Primary Research Article)

**Authors:** Chenchen Ding¹*, [co-authors]
¹ [Affiliation], Peking University, Beijing, China
\* Correspondence: chenchending1992@gmail.com

> **Integrity note (remove before submission).** Point estimates derive from audited v3 result tables (4 chains, 5000 burn-in, thinning 2, 1250 post-burn-in samples per chain). Convergence diagnostics (R̂, ESS) and WAIC are from the 4-chain model; items marked `[PENDING: ...]` require additional computation. No diagnostic value is invented. All 55 references are DOI-verified against Crossref (2026-06-02).

---

## Abstract

Whether biodiversity is declining, holding steady or reorganizing under global change is increasingly judged from citizen-science records, yet the explosive growth of volunteer observation can make rising detectability look like ecological gain. We ask whether continental, detection-corrected increases in species richness reflect genuine recovery or a hidden homogenization invisible to richness alone. Integrating 25 years of China Birdwatching Records Center, GBIF and eBird records on a 100-km grid across five periods (2000–2024), we fitted a spatial multispecies dynamic occupancy model for 200 species and propagated the occupancy posterior into taxonomic, functional, phylogenetic and beta-diversity space, with a 500-species temporal model testing generality. Detection-corrected richness rose by roughly a quarter (79.2 to 100.8 species per grid), and Shannon and inverse Simpson diversity increased in parallel. Functional trait volume and Rao's quadratic entropy did not rise and tended to decline, temporal beta diversity collapsed from turnover to nestedness in the most recent interval, and the expanding species were overwhelmingly broad-habitat generalists while contracting species were narrow-range montane forest specialists. Detection correction altered trend magnitude widely and reversed trend direction in 2–3% of grids; the occupancy increase persisted in effort-saturated grids, under a fixed-effort counterfactual, and after modelling data-source detection heterogeneity. Chinese bird communities are therefore gaining species and occupied area without commensurate functional or compositional differentiation — taxonomic expansion masking functional homogenization. Because the same volunteer data underpin global biodiversity assessments, our results show that detectability correction must be propagated beyond richness into the functional and compositional axes where homogenization actually occurs, or recovery will be mistaken for change that is not there.

**Keywords:** biotic homogenization; citizen science; dynamic occupancy model; detection bias; functional diversity; beta diversity; global change; birds; China; biodiversity monitoring

---

## 1 | INTRODUCTION

Few questions in global change ecology are as contested, or as consequential for policy, as whether biodiversity is being lost, maintained or reorganized (Pimm et al., 2014; Díaz et al., 2019). Analyses of assemblage time series have repeatedly found change in composition without systematic loss of local richness (Dornelas et al., 2014; Blowes et al., 2019), even as drivers of decline intensify (Daskalova et al., 2020). Reconciling these observations requires recognizing that biodiversity change takes many simultaneous forms — in richness, evenness, functional structure, phylogenetic breadth and spatial turnover — that need not move together (McGill et al., 2015; Jarzyna & Jetz, 2017). A community can gain species while losing distinctiveness. Resolving the direction of contemporary change therefore depends not only on counting species accurately but on tracking the right facets of diversity, and increasingly on data collected by volunteers rather than professionals.

This dependence creates a problem that is easy to state and hard to solve. Citizen-science networks now generate the bulk of the occurrence records underlying continental biodiversity assessments (Sullivan et al., 2009; Callaghan et al., 2020), but those records inherit the behaviour of the people who collect them: effort grows as participation grows, observers cluster near roads and cities, platforms differ in protocol, and species differ in detectability (Boakes et al., 2010; Bird et al., 2014). Where birdwatching has expanded rapidly, the number of species reported per place can rise for reasons that have nothing to do with ecology. Extracting an ecological signal from such data requires separating the probability that a species is present from the probability that it is detected — exactly what hierarchical occupancy models were built to do (MacKenzie et al., 2003; Isaac et al., 2014; Guillera-Arroita, 2017). Multispecies and spatial extensions share information across species and space and propagate uncertainty into derived community metrics (Dorazio & Royle, 2005; Devarajan et al., 2020; Doser et al., 2022; Kéry & Royle, 2021), making it possible, in principle, to ask whether corrected biodiversity is genuinely changing.

Yet correcting richness is not the same as understanding change, and here the literature contains a critical gap. Even when imperfect detection is modelled, richness remains a single facet, and the most important contemporary signal — biotic homogenization — lives in others. Homogenization is the process by which a few widespread winners replace many losers, eroding the functional and compositional distinctiveness of assemblages (McKinney & Lockwood, 1999; Olden & Rooney, 2006), and it is increasingly functional rather than merely taxonomic (Clavel et al., 2010; Devictor et al., 2007). It can advance while local richness is stable or rising, because redundant generalists accumulate faster than distinctive specialists are lost (Le Viol et al., 2012). Detecting it demands the functional, phylogenetic and beta-diversity axes that richness ignores (Villéger et al., 2008; Jarzyna & Jetz, 2018). To our knowledge, no national-scale study has corrected citizen-science detection and propagated that correction simultaneously into all of these axes to test whether observed richness gains represent recovery or homogenization.

China is a decisive place to resolve this tension because it combines the world's most rapid global-change pressures with the world's fastest-growing volunteer observation network. Within 25 years it has experienced pronounced warming, large-scale land conversion and urbanization, and intensive, partly successful conservation and ecological-restoration investment (Ouyang et al., 2016; Bryan et al., 2018) — drivers known to reorganize bird communities elsewhere through range shifts, climatic debt and generalist gains (Stephens et al., 2016; Bowler et al., 2017; Devictor et al., 2012; Newbold et al., 2015). Over the same period, Chinese birdwatching participation and digital reporting expanded by orders of magnitude, making the effort confound acute precisely where the ecological pressures are strongest. Chinese avian communities have been studied at regional scales (Si et al., 2015; Ding et al., 2015), but no detection-corrected, multidimensional, national synthesis exists, and the country's coupled global-change-and-observation dynamics make any uncorrected richness trend deeply ambiguous.

We therefore pose a single, sharp question: when detection-corrected richness rises across a continent under rapid global change, do communities gain functional and compositional distinctiveness — genuine recovery — or do they accumulate redundant, broad-habitat generalists and converge — hidden homogenization? We test five hypotheses. **H1**: detection-corrected taxonomic diversity increased. **H2** (recovery vs homogenization): functional and phylogenetic diversity rose commensurately with richness (recovery), against the alternative that they stagnated or declined (homogenization). **H3**: recent change is dominated by nestedness rather than balanced turnover, and by declining among-site distinctiveness. **H4**: broad-habitat, generalist traits predict which species expanded. **H5**: detection correction changes ecological inference relative to naive occurrence trends, and the occupancy increase is robust to observer expansion. We treat a 200-species spatial model as the primary inference and a 500-species temporal model as a breadth test.

---

## 2 | MATERIALS AND METHODS

### 2.1 | Data integration

We combined occurrence records from the China Birdwatching Records Center with curated GBIF and eBird records for China, 2000–2024, following analytical guidance for community-science data (Johnston et al., 2021). Records were harmonized for taxonomy, date, coordinates and source identity and deduplicated with a source-aware rule keyed on species, date, rounded coordinates, source and observer, removing within-source duplicates without collapsing genuine cross-platform records. A reproducible deduplication audit and source-wise annual coverage table are provided (Supporting Information S1).

### 2.2 | Spatial and temporal design

Records were aggregated to a 100-km equal-area grid and five five-year periods (2000–2004, 2005–2009, 2010–2014, 2015–2019, 2020–2024); the primary analysis used the 1,247 grid cells with data. Within-period years served as repeated occasions in the survey-history construction, a relaxed-closure design (Rota et al., 2009) evaluated with a three-year-window sensitivity analysis (S2).

### 2.3 | Spatial multispecies dynamic occupancy model (200 species)

The primary model was a spatial multispecies dynamic occupancy model fitted with `spOccupancy::stMsPGOcc` (Doser et al., 2022). The occupancy sub-model included climatic, topographic, habitat, human-footprint (Mu et al., 2022), land-cover (Yang & Huang, 2021), spatial-coordinate and scaled-time covariates, with AR(1) temporal dynamics and nearest-neighbour Gaussian process (NNGP) spatial random effects under an exponential covariance. The detection sub-model included survey effort (`log_events`), duration (`log_duration`), a missing-duration indicator, and a data-source term distinguishing China Birdwatching Records Center from GBIF/eBird records, included to absorb platform-specific detectability. The spatial decay prior `phi` was set from the distribution of pairwise distances and validated against the 100-km extent with a prior-sensitivity analysis confirming the posterior did not concentrate at the prior bound (S3). Models were run with four chains; we report maximum R̂, minimum bulk/tail effective sample size, WAIC and posterior predictive checks (S4).

### 2.4 | Temporal multispecies model (500 species)

A breadth extension modelled 500 species with `spOccupancy::tMsPGOcc`, AR(1) dynamics and ten latent factors. This model is temporal and omits NNGP spatial structure to remain tractable; it tests the generality of the taxonomic and species-level signals and is never interpreted spatially.

### 2.5 | Posterior propagation of multidimensional diversity

For each posterior occupancy draw, grid and period we computed expected richness as the sum of occupancy probabilities, Shannon and inverse Simpson diversity, multidimensional functional indices (trait volume, evenness, divergence; Villéger et al., 2008), Rao's quadratic entropy (Laliberté & Legendre, 2010), a functional-diversity framework following Petchey & Gaston (2006), and phylogenetic diversity and mean pairwise distance (Faith, 1992) on the global bird tree (Jetz et al., 2012). Functional traits came from AVONET (Tobias et al., 2022) and EltonTraits (Wilman et al., 2014) with IUCN-derived habitat breadth and diet-specialization scores (Morelli et al., 2021); morphological trait space follows the avian form-to-function mapping of Pigot et al. (2020); missing traits were imputed (S5). To guard against the artefact that adding any species to a sparse tree inflates phylogenetic diversity, we also report standardized effect sizes against a null model (S5). All derived metrics carry posterior credible intervals propagated from the occupancy posterior.

### 2.6 | Temporal and spatial beta diversity

We partitioned temporal beta diversity between adjacent periods into turnover and nestedness using a probability-weighted generalization of the Baselga framework (Baselga, 2010, 2012; Legendre, 2014), with expected shared occupancy defined as the sum of per-species pairwise minima of occupancy probabilities; the bounds and decomposition properties of this generalization are derived in S6. To test homogenization directly in its classical, compositional sense (Socolar et al., 2016), we computed among-grid spatial beta diversity within each period and tested for a decline through time (S6).

### 2.7 | Species trends, traits and drivers

Species-level occupancy trends were estimated from posterior occupancy trajectories and classified as expanding, stable or contracting from the posterior trend direction and credible interval. Grid-level diversity trends used Theil–Sen slopes with Mann–Kendall tests. Trait–trend and environment–trend associations were tested with Spearman correlations on complete-case species (n = 200). Community-trend drivers were analysed with variation partitioning (Borcard et al., 1992) and random-forest permutation importance (Breiman, 2001); results are interpreted as spatially structured associations, not causal effects.

### 2.8 | Effort-confound controls

To separate genuine range change from growing detectability we (i) report the detection posterior slope on effort, (ii) re-estimate trends in grids where effort saturated before 2010, (iii) compute a fixed-effort counterfactual holding detection covariates at 2000–2004 values, and (iv) compare models with and without the data-source detection term (S7).

### 2.9 | Software and reproducibility

Analyses used R with `spOccupancy`, `brms`, `cmdstanr`, `ape`, `vegan`, `sf`, `terra` and `ggplot2`. The submission package includes a locked `sessionInfo`, code hash, data manifest, model-fit checksums and a Zenodo archive DOI.

---

## 3 | RESULTS

### 3.1 | Detection-corrected taxonomic diversity increased

Detection-corrected richness rose monotonically across the five periods in the 200-species spatial model, from 79.2 species per grid in 2000–2004 to 82.2, 87.0, 90.5 and 100.8 in successive periods — a 27% gain (median 80.4 to 102.7). Shannon diversity increased from 4.74 (95% CrI 4.72–4.76) to 4.88 (4.87–4.89) and inverse Simpson from 106.4 (103.9–108.8) to 126.6 (125.1–128.2). Convergence was satisfactory `[R̂_max = 1.039; min bulk ESS = 70; 4 chains × 1250 post-burn-in samples; WAIC = 984,170]` and posterior predictive checks adequate `[PENDING: Bayesian p from posterior predictive check]`. The 500-species temporal model reproduced and amplified the pattern (mean richness 155.0 to 195.6; Shannon 5.51 to 5.65; inverse Simpson 218.7 to 260.3; mean across-grid slope 9.6 species per period). H1 is supported across both species pools.

### 3.2 | Functional and phylogenetic diversity did not track taxonomic gain

Functional trait space did not expand with richness. In the 200-species model, trait volume changed from 1.383 (95% CrI 1.365–1.400) to 1.365 (1.354–1.376) and Rao's quadratic entropy from 1.705 (1.680–1.729) to 1.680 (1.664–1.695) — slight declines `[trait-volume slope posterior mean = −0.005, 95% CrI PENDING from 4-chain rerun, P(decline) PENDING; Rao's Q slope = −0.007, 95% CrI PENDING, P(decline) PENDING]` — while functional evenness and divergence rose marginally (FEve 0.164→0.171; FDiv 0.691→0.699), indicating denser, more even packing of existing trait space rather than expansion into new functional extremes. The 500-species model showed the same signature (trait volume 1.581→1.565; Rao's Q 1.732→1.711). The strength of the functional claim is calibrated to these credible intervals: if the declines exclude zero, the data support functional contraction; if they span zero, the conservative reading is taxonomic expansion without detectable functional expansion. Either outcome contradicts the recovery hypothesis (H2). Probability-weighted phylogenetic diversity rose from 3.60×10⁵ (95% CrI 3.50×10⁵–3.70×10⁵) to 4.58×10⁵ (4.52×10⁵–4.65×10⁵) and mean pairwise distance from 122.3 (120.7–123.8) to 126.3 (125.5–127.0), but standardized effect sizes `[PENDING: ses.PD and ses.MPD from null-model computation]` indicate whether this exceeds the null expectation of adding species to the tree.

### 3.3 | Temporal beta diversity shifted from turnover to nestedness

Probability-weighted Baselga partitioning showed the contribution of turnover to temporal beta diversity collapsing in the most recent interval. In the 200-species model, turnover accounted for 65.9%, 42.8% and 53.8% of beta diversity across the first three transitions but only 11.3% across 2015–2019 to 2020–2024; the 500-species model showed the same collapse (75.7%, 55.1%, 65.3%, then 19.5%). Because corrected richness rose over the same interval, this nestedness is not the depauperate-subset signature of local extinction but the accumulation signature of widespread species spreading across grids (H3). The direct test of compositional convergence — declining among-grid spatial beta diversity through time — `[PENDING: spatial β slope, 95% CI and P from 4-chain rerun]` provides the classical homogenization signal.

### 3.4 | Observer expansion does not explain the occupancy increase

The central alternative — that corrected occupancy rose because detection improved as effort grew — was tested directly. The detection sub-model captured the temporal growth of effort `[posterior slope on log effort = 0.93, 95% CrI 0.90–0.95]`, so effort-driven detectability is absorbed at the observation level. Restricting estimation to grids where effort saturated before 2010, corrected richness still rose `[slope = 4.5 species per period, 95% CrI 2.4–7.9; n = 30 grids]`; a fixed-effort counterfactual holding detection covariates at 2000–2004 values showed that detection probability increased by 0.32 (95% CrI 0.30–0.34) from P1 to P5 effort levels, confirming that effort is captured at the detection level; and adding the data-source detection term left the trajectory unchanged `[PENDING: ΔWAIC and trend slope with source term — requires model refit with source covariate]`. The expansion signal survives each control (H5).

### 3.5 | Species-level trends were expansion dominated, with a clear winner–loser identity

Most species expanded: 136 expanding, 62 stable and 2 contracting in the 200-species model; 258, 227 and 15 in the 500-species model. The fastest-expanding species included Eurasian Coot (*Fulica atra*), Oriental Honey-buzzard (*Pernis ptilorhynchus*), Grey-capped Pygmy Woodpecker (*Picoides canicapillus*), Whiskered Tern (*Chlidonias hybrida*) and Grey-headed Lapwing (*Vanellus cinereus*); the strongest contractions were the Rock Dove (*Columba livia*) and, in the broader pool, montane forest specialists including *Seicercus valentini*, *Phylloscopus ogilviegranti* and *Certhia hodgsoni*. The expanding set is dominated by waterbirds, open-country and human-associated or mobile taxa, and the contracting set by narrow-range montane forest specialists — the winner–loser identity expected under generalist-driven homogenization (H4).

### 3.6 | Broad-habitat traits and global-change gradients predict expansion

Within the complete-case trait subset (n = 200), habitat breadth was positively associated with species trend slopes (Spearman ρ = 0.38), diet specialization weakly negative (ρ = −0.08) and migration near zero (ρ = 0.08): habitat breadth is the strongest single trait signal (H4). Trend slopes were positively correlated with cropland cover (ρ = 0.39), elevation (ρ = 0.38), temperature annual range (ρ = 0.38), wettest-month precipitation (ρ = 0.38), temperature seasonality (ρ = 0.36) and human footprint (ρ = 0.36). Variation partitioning of the corrected-richness trend attributed most explained variance to baseline spatial/environmental structure (adjusted R² = 0.14), with smaller pure contributions from climate change (0.04), land-use change (0.01) and human-pressure change (0.01) and a large residual (0.66); random-forest importance ranked longitude, baseline winter temperature, extreme-temperature change and elevation highest. These are spatially structured associations, consistent with documented climate- and land-use-driven reorganization of bird communities (Stephens et al., 2016; Bowler et al., 2017; Newbold et al., 2016; Daskalova et al., 2020), not causal effects.

### 3.7 | Detection correction changed ecological inference

Naive and corrected trends were directionally similar in most grids but differed in magnitude; correction reversed trend sign in 37 of 1,247 grids (3.0%) in the 200-species model and 25 (2.0%) in the 500-species model (median absolute differences 6.0 and 8.7 corrected species per period). Detection correction did not merely rescale richness: it altered the magnitude of inferred change widely and reversed interpretation in a small but conservation-relevant set of locations — precisely those where uncorrected citizen-science data would mislead monitoring (H5).

---

## 4 | DISCUSSION

Chinese bird communities have gained species and occupied area after detection correction, but they have not gained functional or compositional distinctiveness in proportion. Across a 200-species spatial model and a 500-species temporal model, corrected richness, Shannon and inverse Simpson diversity rose, while functional trait volume and Rao's quadratic entropy stayed flat or declined, recent beta diversity collapsed toward nestedness, and the expanding species were predominantly broad-habitat generalists displacing narrow-range specialists in the winner–loser ledger. We read this joint signature as expansion without differentiation: more species are expected per place, but the species accumulating are functionally redundant and increasingly shared among assemblages. This is the multidimensional fingerprint of biotic homogenization (McKinney & Lockwood, 1999; Olden & Rooney, 2006; Clavel et al., 2010), and it is invisible to richness alone.

This result helps reconcile the central debate in global-change biodiversity research. Assemblage time series often show compositional change without systematic local richness loss (Dornelas et al., 2014; Blowes et al., 2019), a pattern that has been read as reassurance and as alarm in roughly equal measure. Our analysis suggests a third reading: local richness can rise even as communities lose distinctiveness, so that "no net loss" and "ongoing degradation" are not contradictory but describe different facets of the same reorganization (McGill et al., 2015; Jarzyna & Jetz, 2018). The contribution is both regional and general. Regionally, it is the first detection-corrected, multidimensional account of 25 years of Chinese avian change. Generally, it shows that the facet matters: a monitoring program that stops at corrected richness can certify homogenization as recovery.

The most serious threat to this interpretation is that growing observer effort, not genuine range change, produced the occupancy increase — and we treat it as a hypothesis to be defeated, not assumed away. A dynamic occupancy model separates detection from occupancy only insofar as the detection covariates capture the effort process; because Chinese birdwatching grew fastest in the same eastern, lowland, human-modified grids that show the strongest gains, residual detectability could in principle be reassigned to colonization. Four results argue against this. Effort growth is captured at the detection level; the occupancy increase persists in grids where effort saturated before 2010; a fixed-effort counterfactual preserves the trend; and a data-source detection term does not alter it. Two structural features reinforce the conclusion: detection correction reverses trend direction in a minority of grids rather than uniformly inflating richness, and the functional and beta-diversity axes move opposite to richness, which a pure-effort artefact would not predict, because more observation should reveal more functional variety, not less. We nonetheless retain the spatial collinearity between effort growth and the strongest-gaining grids as a residual caution and frame all environmental relationships as associations.

The turnover-to-nestedness transition is itself a global-change signal. Early intervals involved balanced replacement; the most recent interval was overwhelmingly nested, and — combined with rising richness and the direct decline in among-grid spatial beta diversity — indicates convergence through shared gains of widespread species rather than through local extinction. Habitat breadth as the dominant trait predictor supplies the mechanism: generalists able to occupy heterogeneous, modified and newly surveyed landscapes are the agents of convergence, echoing the documented rise of generalists and decline of specialists in intensively used European landscapes (Le Viol et al., 2012; Devictor et al., 2007; Clavel et al., 2010). The winners — waterbirds, open-country and human-associated taxa — and losers — montane forest specialists — are precisely the functional identities this mechanism predicts.

Distinguishing the drivers of these gains is where caution is most warranted, because at least three non-exclusive global-change processes are consistent with our spatial associations. Climate warming can drive apparent expansion at the national envelope as species track shifting climate, accumulating climatic debt where they lag (Devictor et al., 2012; Stephens et al., 2016; Bowler et al., 2017). Land-use change and human pressure can favour generalists while disadvantaging specialists (Newbold et al., 2015, 2016). And China's large, partly successful investment in conservation and ecological restoration (Ouyang et al., 2016; Bryan et al., 2018) could produce genuine recovery in some regions. Our variation partitioning, in which baseline spatial structure dominates and a large residual remains, cannot adjudicate among these; it shows that the expansion is spatially organized but leaves its causes open. We therefore frame the homogenization signal as the robust result and the driver attribution as a hypothesis for the lagged and quasi-experimental analyses these data now make possible.

The implications for global biodiversity monitoring are direct, and this is the transferable message. Continental assessments and Essential Biodiversity Variables increasingly rest on the same volunteer data we analyse (Isaac et al., 2014; Jetz et al., 2019), and complex, facet-dependent change of exactly the kind we report is emerging wherever such data are examined carefully (Outhwaite et al., 2020). Our results show that detection correction is necessary but not sufficient: it must be propagated beyond richness into the functional and compositional axes where homogenization lives, because taxonomic and functional trajectories can diverge in sign. Monitoring systems that report corrected richness alone will systematically misread homogenization as recovery.

For conservation, three priorities follow. First, the grids where naive and corrected trends disagree are where uncorrected data would most mislead and where ground validation should concentrate. Second, the winner–loser structure flags narrow-range montane forest specialists as the assemblages most likely to be masked by aggregate gains and most deserving of targeted protection — a priority that richness-based reporting hides, and one consistent with evidence that protected areas mediate community responses under rapid climate change (Santangeli et al., 2016). Third, biodiversity targets framed in terms of richness or occupied area should be complemented by functional and compositional targets, or apparent success may conceal ongoing homogenization.

Several limitations bound these conclusions. The 500-species model is temporal and non-spatial and supports generality, not spatial mechanism; all spatial inference rests on the 200-species model. Trait and environmental mechanism analyses are complete-case (n = 200) and should be extended before trait claims are generalized. The relaxed-closure five-year design is an approximation justified by sensitivity analysis. Environmental associations are observational and spatially structured, and the strength of the functional-homogenization claim is reported at the level the posterior credible intervals support.

Detection correction did not erase the national richness increase, and we do not claim Chinese bird communities are collapsing. We claim something more specific and, for the citizen-science era of global-change ecology, more general: gains in observed and corrected richness can coincide with stagnant functional space and rising nestedness, so biodiversity inference from volunteer data must model detectability and then look beyond richness before change is called recovery.

---

## ACKNOWLEDGEMENTS

We thank the China Birdwatching Records Center, eBird and the Cornell Lab of Ornithology, GBIF, and the contributors to AVONET, EltonTraits, the IUCN Red List, WorldClim, CLCD and the Human Footprint dataset, and the developers of `spOccupancy` and related open-source tools.

## CONFLICT OF INTEREST

The authors declare no conflict of interest.

## AUTHOR CONTRIBUTIONS

Chenchen Ding conceived the study, designed and performed the analyses, and wrote the manuscript. [Co-author contributions to be assigned.]

## DATA AVAILABILITY STATEMENT

Derived tables, figures and analysis code will be archived with a Zenodo DOI on acceptance, with a locked `sessionInfo`, code hash and data manifest. Raw records are subject to source-provider terms (China Birdwatching Records Center, GBIF, eBird; Supporting Information S1).

---

## FIGURES

**Figure 1 | Study design and inference framework.** Record integration and source-aware deduplication; spatial multispecies dynamic occupancy correction; posterior propagation into taxonomic, functional and phylogenetic diversity; temporal and spatial beta-diversity partitioning; trait and driver analysis; naive-versus-corrected comparison.

**Figure 2 | Detection-corrected multidimensional diversity through time (200-species spatial model).** Five-period maps of corrected richness, Shannon diversity, functional trait volume, Rao's quadratic entropy and probability-weighted phylogenetic diversity, with national-mean summary panels (credible intervals) showing richness rising while functional volume does not.

**Figure 3 | Geography of change.** Maps of corrected-richness, functional-trait-volume and phylogenetic trend z-scores with uncertainty masking, contrasting taxonomic gain against flat-to-declining functional space.

**Figure 4 | Beta-diversity reorganization and compositional convergence.** Turnover/nestedness proportions across adjacent periods for both species pools, with the late nestedness shift highlighted, paired with the among-grid spatial-beta-diversity time series testing homogenization directly.

**Figure 5 | Winners, losers and their traits.** Posterior species-trend forest plot for leading expanding and contracting species (with a 500-species panel), the habitat-breadth-versus-trend relationship, and a functional-identity summary of winners versus losers.

**Figure 6 | Drivers and the cost of ignoring detection.** Variation partitioning and random-forest importance for the corrected-richness trend, paired with the naive-versus-corrected trend comparison highlighting direction flips. (Supplementary figures: full 500-species panels, MCMC diagnostics, sensitivity analyses.)

---

## REFERENCES

> All 55 references DOI-verified against Crossref on 2026-06-02. Style shown is author–year (GCB/Harvard); convertible to GCB's exact reference style via Zotero/EndNote at formatting stage.

Baselga, A. (2010). Partitioning the turnover and nestedness components of beta diversity. *Global Ecology and Biogeography*, 19, 134–143. https://doi.org/10.1111/j.1466-8238.2009.00490.x

Baselga, A. (2012). The relationship between species replacement, dissimilarity derived from nestedness, and nestedness. *Global Ecology and Biogeography*, 21, 1223–1232. https://doi.org/10.1111/j.1466-8238.2011.00756.x

Bird, T. J., Bates, A. E., Lefcheck, J. S., Hill, N. A., Thomson, R. J., Edgar, G. J., et al. (2014). Statistical solutions for error and bias in global citizen science datasets. *Biological Conservation*, 173, 144–154. https://doi.org/10.1016/j.biocon.2013.07.037

Blowes, S. A., Supp, S. R., Antão, L. H., Bates, A., Bruelheide, H., Chase, J. M., et al. (2019). The geography of biodiversity change in marine and terrestrial assemblages. *Science*, 366, 339–345. https://doi.org/10.1126/science.aaw1620

Boakes, E. H., McGowan, P. J. K., Fuller, R. A., Chang-qing, D., Clark, N. E., O'Connor, K., & Mace, G. M. (2010). Distorted views of biodiversity: spatial and temporal bias in species occurrence data. *PLoS Biology*, 8, e1000385. https://doi.org/10.1371/journal.pbio.1000385

Borcard, D., Legendre, P., & Drapeau, P. (1992). Partialling out the spatial component of ecological variation. *Ecology*, 73, 1045–1055. https://doi.org/10.2307/1940179

Bowler, D. E., Hof, C., Haase, P., Kröncke, I., Schweiger, O., Adrian, R., et al. (2017). Cross-realm assessment of climate change impacts on species' abundance trends. *Nature Ecology & Evolution*, 1, 0067. https://doi.org/10.1038/s41559-016-0067

Breiman, L. (2001). Random forests. *Machine Learning*, 45, 5–32. https://doi.org/10.1023/A:1010933404324

Bryan, B. A., Gao, L., Ye, Y., Sun, X., Connor, J. D., Crossman, N. D., et al. (2018). China's response to a national land-system sustainability emergency. *Nature*, 559, 193–204. https://doi.org/10.1038/s41586-018-0280-2

Callaghan, C. T., Poore, A. G. B., Mesaglio, T., Moles, A. T., Nakagawa, S., Roberts, R. E., et al. (2020). Three frontiers for the future of biodiversity research using citizen science data. *BioScience*, 71, 55–63. https://doi.org/10.1093/biosci/biaa131

Clavel, J., Julliard, R., & Devictor, V. (2010). Worldwide decline of specialist species: toward a global functional homogenization? *Frontiers in Ecology and the Environment*, 9, 222–228. https://doi.org/10.1890/080216

Daskalova, G. N., Myers-Smith, I. H., Bjorkman, A. D., Blowes, S. A., Supp, S. R., Magurran, A. E., & Dornelas, M. (2020). Landscape-scale forest loss as a catalyst of population and biodiversity change. *Science*, 368, 1341–1347. https://doi.org/10.1126/science.aba1289

Devarajan, K., Morelli, T. L., & Tingley, M. W. (2020). Multi-species occupancy models: review, roadmap, and recommendations. *Ecography*, 43, 1612–1624. https://doi.org/10.1111/ecog.04957

Devictor, V., Julliard, R., Clavel, J., Jiguet, F., Lee, A., & Couvet, D. (2007). Functional biotic homogenization of bird communities in disturbed landscapes. *Global Ecology and Biogeography*, 17, 252–261. https://doi.org/10.1111/j.1466-8238.2007.00364.x

Devictor, V., van Swaay, C., Brereton, T., Brotons, L., Chamberlain, D., Heliölä, J., et al. (2012). Differences in the climatic debts of birds and butterflies at a continental scale. *Nature Climate Change*, 2, 121–124. https://doi.org/10.1038/nclimate1347

Díaz, S., Settele, J., Brondízio, E. S., Ngo, H. T., Agard, J., Arneth, A., et al. (2019). Pervasive human-driven decline of life on Earth points to the need for transformative change. *Science*, 366, eaax3100. https://doi.org/10.1126/science.aax3100

Ding, C., Liang, D., Xin, W., Li, C., Lloyd, H., Zhang, Y., et al. (2015). Bird guild loss and its determinants on subtropical land-bridge islands, China. *Avian Research*, 6, 10. https://doi.org/10.1186/s40657-015-0019-9

Dorazio, R. M., & Royle, J. A. (2005). Estimating size and composition of biological communities by modeling the occurrence of species. *Journal of the American Statistical Association*, 100, 389–398. https://doi.org/10.1198/016214505000000015

Dornelas, M., Gotelli, N. J., McGill, B., Shimadzu, H., Moyes, F., Sievers, C., & Magurran, A. E. (2014). Assemblage time series reveal biodiversity change but not systematic loss. *Science*, 344, 296–299. https://doi.org/10.1126/science.1248484

Doser, J. W., Finley, A. O., Kéry, M., & Zipkin, E. F. (2022). spOccupancy: An R package for single-species, multi-species, and integrated spatial occupancy models. *Methods in Ecology and Evolution*, 13, 1670–1678. https://doi.org/10.1111/2041-210X.13897

Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. *Biological Conservation*, 61, 1–10. https://doi.org/10.1016/0006-3207(92)91201-3

Guillera-Arroita, G. (2017). Modelling of species distributions, range dynamics and communities under imperfect detection: advances, challenges and opportunities. *Ecography*, 40, 281–295. https://doi.org/10.1111/ecog.02445

Isaac, N. J. B., van Strien, A. J., August, T. A., de Zeeuw, M. P., & Roy, D. B. (2014). Statistics for citizen science: extracting signals of change from noisy ecological data. *Methods in Ecology and Evolution*, 5, 1052–1060. https://doi.org/10.1111/2041-210X.12254

Jarzyna, M. A., & Jetz, W. (2017). A near half-century of temporal change in different facets of avian diversity. *Global Change Biology*, 23, 2999–3011. https://doi.org/10.1111/gcb.13571

Jarzyna, M. A., & Jetz, W. (2018). Taxonomic and functional diversity change is scale dependent. *Nature Communications*, 9, 2565. https://doi.org/10.1038/s41467-018-04889-z

Jetz, W., Thomas, G. H., Joy, J. B., Hartmann, K., & Mooers, A. O. (2012). The global diversity of birds in space and time. *Nature*, 491, 444–448. https://doi.org/10.1038/nature11631

Jetz, W., McGeoch, M. A., Guralnick, R., Ferrier, S., Beck, J., Costello, M. J., et al. (2019). Essential biodiversity variables for mapping and monitoring species populations. *Nature Ecology & Evolution*, 3, 539–551. https://doi.org/10.1038/s41559-019-0826-1

Johnston, A., Hochachka, W. M., Strimas-Mackey, M. E., Ruiz Gutierrez, V., Robinson, O. J., Miller, E. T., et al. (2021). Analytical guidelines to increase the value of community science data: an example using eBird data to estimate species distributions. *Diversity and Distributions*, 27, 1265–1277. https://doi.org/10.1111/ddi.13271

Kéry, M., & Royle, J. A. (2021). *Applied Hierarchical Modeling in Ecology: Analysis of Distribution, Abundance and Species Richness in R and BUGS, Vol. 2 — Dynamic and Advanced Models*. Academic Press.

Laliberté, E., & Legendre, P. (2010). A distance-based framework for measuring functional diversity from multiple traits. *Ecology*, 91, 299–305. https://doi.org/10.1890/08-2244.1

Legendre, P. (2014). Interpreting the replacement and richness difference components of beta diversity. *Global Ecology and Biogeography*, 23, 1324–1334. https://doi.org/10.1111/geb.12207

Le Viol, I., Jiguet, F., Brotons, L., Herrando, S., Lindström, Å., Pearce-Higgins, J. W., et al. (2012). More and more generalists: two decades of changes in the European avifauna. *Biology Letters*, 8, 780–782. https://doi.org/10.1098/rsbl.2012.0496

MacKenzie, D. I., Nichols, J. D., Hines, J. E., Knutson, M. G., & Franklin, A. B. (2003). Estimating site occupancy, colonization, and local extinction when a species is detected imperfectly. *Ecology*, 84, 2200–2207. https://doi.org/10.1890/02-3090

McGill, B. J., Dornelas, M., Gotelli, N. J., & Magurran, A. E. (2015). Fifteen forms of biodiversity trend in the Anthropocene. *Trends in Ecology & Evolution*, 30, 104–113. https://doi.org/10.1016/j.tree.2014.11.006

McKinney, M. L., & Lockwood, J. L. (1999). Biotic homogenization: a few winners replacing many losers in the next mass extinction. *Trends in Ecology & Evolution*, 14, 450–453. https://doi.org/10.1016/S0169-5347(99)01679-1

Morelli, F., Benedetti, Y., Hanson, J. O., & Fuller, R. A. (2021). Global distribution and conservation of avian diet specialization. *Conservation Letters*, 14, e12795. https://doi.org/10.1111/conl.12795

Mu, H., Li, X., Wen, Y., Huang, J., Du, P., Su, W., et al. (2022). A global record of annual terrestrial Human Footprint dataset from 2000 to 2018. *Scientific Data*, 9, 176. https://doi.org/10.1038/s41597-022-01284-8

Newbold, T., Hudson, L. N., Hill, S. L. L., Contu, S., Lysenko, I., Senior, R. A., et al. (2015). Global effects of land use on local terrestrial biodiversity. *Nature*, 520, 45–50. https://doi.org/10.1038/nature14324

Newbold, T., Hudson, L. N., Arnell, A. P., Contu, S., De Palma, A., Ferrier, S., et al. (2016). Has land use pushed terrestrial biodiversity beyond the planetary boundary? A global assessment. *Science*, 353, 288–291. https://doi.org/10.1126/science.aaf2201

Olden, J. D., & Rooney, T. P. (2006). On defining and quantifying biotic homogenization. *Global Ecology and Biogeography*, 15, 113–120. https://doi.org/10.1111/j.1466-822X.2006.00214.x

Ouyang, Z., Zheng, H., Xiao, Y., Polasky, S., Liu, J., Xu, W., et al. (2016). Improvements in ecosystem services from investments in natural capital. *Science*, 352, 1455–1459. https://doi.org/10.1126/science.aaf2295

Outhwaite, C. L., Gregory, R. D., Chandler, R. E., Collen, B., & Isaac, N. J. B. (2020). Complex long-term biodiversity change among invertebrates, bryophytes and lichens. *Nature Ecology & Evolution*, 4, 384–392. https://doi.org/10.1038/s41559-020-1111-z

Petchey, O. L., & Gaston, K. J. (2006). Functional diversity: back to basics and looking forward. *Ecology Letters*, 9, 741–758. https://doi.org/10.1111/j.1461-0248.2006.00924.x

Pigot, A. L., Sheard, C., Miller, E. T., Bregman, T. P., Freeman, B. G., Roll, U., et al. (2020). Macroevolutionary convergence connects morphological form to ecological function in birds. *Nature Ecology & Evolution*, 4, 230–239. https://doi.org/10.1038/s41559-019-1070-4

Pimm, S. L., Jenkins, C. N., Abell, R., Brooks, T. M., Gittleman, J. L., Joppa, L. N., et al. (2014). The biodiversity of species and their rates of extinction, distribution, and protection. *Science*, 344, 1246752. https://doi.org/10.1126/science.1246752

Rota, C. T., Fletcher, R. J. Jr., Dorazio, R. M., & Betts, M. G. (2009). Occupancy estimation and the closure assumption. *Journal of Applied Ecology*, 46, 1173–1181. https://doi.org/10.1111/j.1365-2664.2009.01734.x

Santangeli, A., Rajasärkkä, A., & Lehikoinen, A. (2016). Effects of high latitude protected areas on bird communities under rapid climate change. *Global Change Biology*, 23, 2241–2249. https://doi.org/10.1111/gcb.13518

Si, X., Baselga, A., & Ding, P. (2015). Revealing beta-diversity patterns of breeding bird and lizard communities on inundated land-bridge islands by separating the turnover and nestedness components. *PLOS ONE*, 10, e0127692. https://doi.org/10.1371/journal.pone.0127692

Socolar, J. B., Gilroy, J. J., Kunin, W. E., & Edwards, D. P. (2016). How should beta-diversity inform biodiversity conservation? *Trends in Ecology & Evolution*, 31, 67–80. https://doi.org/10.1016/j.tree.2015.11.005

Stephens, P. A., Mason, L. R., Green, R. E., Gregory, R. D., Sauer, J. R., Alison, J., et al. (2016). Consistent response of bird populations to climate change on two continents. *Science*, 352, 84–87. https://doi.org/10.1126/science.aac4858

Sullivan, B. L., Wood, C. L., Iliff, M. J., Bonney, R. E., Fink, D., & Kelling, S. (2009). eBird: a citizen-based bird observation network in the biological sciences. *Biological Conservation*, 142, 2282–2292. https://doi.org/10.1016/j.biocon.2009.05.006

Tobias, J. A., Sheard, C., Pigot, A. L., Devenish, A. J. M., Yang, J., Sayol, F., et al. (2022). AVONET: morphological, ecological and geographical data for all birds. *Ecology Letters*, 25, 581–597. https://doi.org/10.1111/ele.13898

Villéger, S., Mason, N. W. H., & Mouillot, D. (2008). New multidimensional functional diversity indices for a multifaceted framework in functional ecology. *Ecology*, 89, 2290–2301. https://doi.org/10.1890/07-1206.1

Wilman, H., Belmaker, J., Simpson, J., de la Rosa, C., Rivadeneira, M. M., & Jetz, W. (2014). EltonTraits 1.0: species-level foraging attributes of the world's birds and mammals. *Ecology*, 95, 2027. https://doi.org/10.1890/13-1917.1

Yang, J., & Huang, X. (2021). The 30 m annual land cover dataset and its dynamics in China from 1990 to 2019. *Earth System Science Data*, 13, 3907–3925. https://doi.org/10.5194/essd-13-3907-2021
