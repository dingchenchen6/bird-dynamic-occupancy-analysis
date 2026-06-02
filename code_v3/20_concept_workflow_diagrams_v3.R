#!/usr/bin/env Rscript
## 20_concept_workflow_diagrams_v3.R  —  v3 概念 Venn + 工作流 DAG
suppressPackageStartupMessages({
  library(ggplot2); library(ggforce); library(tibble); library(dplyr)
  library(DiagrammeR); library(DiagrammeRsvg); library(rsvg)
})
CODE_V3 <- Sys.getenv("V3_CODE_DIR", "/home/dingchenchen/bird_dynamic_occupancy_analysis/code_v3")
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
source(file.path(CODE_V3, "utils_mapping.R"))
ensure_v3_dirs()
log_time("20", "concept + workflow diagrams")

# ---- 概念 Venn ----
circles <- tibble(
  dim = c("Taxonomic","Phylogenetic","Functional"),
  x = c(-1.3, 1.3, 0), y = c(-0.7, -0.7, 1.0), r = 2.2,
  fill = c("#0E5A78","#6A4C93","#C9784A"))
indices <- tibble(
  dim = c(rep("Taxonomic",3), rep("Phylogenetic",5), rep("Functional",6)),
  label = c("Richness","Shannon","Inverse Simpson",
            "Faith PD","MPD","MNTD","NRI","NTI",
            "FRic","FEve","FDiv","FDis","Rao Q","CWM/CWSD"),
  x = c(-3.0,-3.0,-3.0, 2.4,2.4,2.4,2.4,2.4, -0.5,0.4,-0.5,0.4,-0.5,0.4),
  y = c(0.4,0.0,-0.4, 0.4,0.1,-0.2,-0.5,-0.8, 2.6,2.35,2.1,1.85,1.6,1.35))
p_concept <- ggplot() +
  ggforce::geom_circle(data=circles, aes(x0=x,y0=y,r=r,fill=fill), alpha=0.18, colour=NA) +
  ggforce::geom_circle(data=circles, aes(x0=x,y0=y,r=r,colour=fill), fill=NA, linewidth=0.6) +
  geom_text(data=circles, aes(x=x, y=y+r-0.18, label=dim, colour=fill),
            size=4.5, fontface="bold") +
  geom_text(data=indices, aes(x=x, y=y, label=label,
            colour=circles$fill[match(dim,circles$dim)]), size=2.8, hjust=0) +
  scale_fill_identity() + scale_colour_identity() +
  annotate("text", x=0, y=-2.3, label="Multidimensional community structure",
           size=5, fontface="bold") +
  coord_fixed(xlim=c(-4.2,4.2), ylim=c(-3.2,4.6)) +
  theme_void(base_family="Arial") +
  theme(plot.title=element_text(face="bold", size=12, hjust=0.5),
        plot.margin=margin(8,8,8,8)) +
  labs(title="Multidimensional diversity of bird communities")
# 顶部驱动
drv <- tibble(xmin=c(-3.7,-1.7,0.3,2.3), xmax=c(-1.9,-0.3,1.7,3.7),
              ymin=3.3, ymax=4.1,
              label=c("Climate","Land use","Topography","Human pressure"))
p_concept <- p_concept +
  geom_rect(data=drv, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            fill="#F2E8D8", colour="#7C7C7C", linewidth=0.4) +
  geom_text(data=drv, aes(x=(xmin+xmax)/2,y=(ymin+ymax)/2,label=label),
            size=3, fontface="bold") +
  annotate("segment", x=0, y=3.25, xend=0, yend=2.05,
           colour="grey40", linewidth=0.4,
           arrow=grid::arrow(length=grid::unit(0.16,"cm"), type="closed")) +
  annotate("text", x=0.25, y=2.85, label="drives via\noccupancy psi",
           size=2.6, colour="grey25", hjust=0)
save_nature(p_concept, "fig_v3_concept_diversity", width_mm=183, height_mm=160)

# ---- 工作流 DAG ----
dot <- "digraph workflow {
  graph [layout=dot, rankdir=LR, fontname=Arial, bgcolor=white, pad=0.4, nodesep=0.4, ranksep=0.55]
  node  [fontname=Arial, shape=box, style=\"rounded,filled\", fontsize=10, penwidth=0.6, color=\"#333333\"]
  edge  [fontname=Arial, fontsize=8, color=\"#666666\", penwidth=0.8, arrowsize=0.6]
  subgraph cluster_data { label=\"Data integration\"; color=\"#0E5A78\"; style=\"rounded,dashed\"
    eBird [label=\"eBird / GBIF\", fillcolor=\"#E6F0F4\"]
    BW [label=\"China Birdwatch\", fillcolor=\"#E6F0F4\"]
    iNat [label=\"iNaturalist\", fillcolor=\"#E6F0F4\"]
    Dedup [label=\"Deduplication\", fillcolor=\"#F0DCC0\"]
  }
  subgraph cluster_grid { label=\"Spatial-temporal frame\"; color=\"#6A4C93\"; style=\"rounded,dashed\"
    Grid [label=\"100 km grid\", fillcolor=\"#E9DEF2\"]
    Blocks [label=\"5-yr periods\\n2000-2024\", fillcolor=\"#E9DEF2\"]
    Cov [label=\"Detection covariates\", fillcolor=\"#E9DEF2\"]
  }
  subgraph cluster_model { label=\"Bayesian model\"; color=\"#C9784A\"; style=\"rounded,dashed\"
    Fit [label=\"stMsPGOcc\\n4 chains, AR1, spatial GP\", fillcolor=\"#F8E0D2\", fontname=\"Arial-Bold\"]
    Diag [label=\"R-hat / ESS / PPC\", fillcolor=\"#F8E0D2\"]
    Psi [label=\"Posterior psi\\n200 draws\", fillcolor=\"#F8E0D2\"]
  }
  subgraph cluster_div { label=\"Multidimensional diversity\"; color=\"#3C8C5A\"; style=\"rounded,dashed\"
    Tax [label=\"Taxonomic\", fillcolor=\"#DCEFE2\"]
    Phy [label=\"Phylogenetic\\nMcTavish tree\", fillcolor=\"#DCEFE2\"]
    Fun [label=\"Functional\\nfundiversity\", fillcolor=\"#DCEFE2\"]
    Beta [label=\"Temporal beta\\nBaselga + homogenisation\", fillcolor=\"#DCEFE2\"]
  }
  subgraph cluster_inf { label=\"Driver inference\"; color=\"#8B2E1E\"; style=\"rounded,dashed\"
    VIF [label=\"VIF + corr\", fillcolor=\"#F5E0E0\"]
    Brms [label=\"brms + cmdstanr\\ngp(lon,lat)\\ninteractions\", fillcolor=\"#F5E0E0\"]
    Imp [label=\"LMG + RF + varpart\", fillcolor=\"#F5E0E0\"]
  }
  subgraph cluster_out { label=\"Outputs\"; color=\"#1F1F1F\"; style=\"rounded,dashed\"
    Maps [label=\"Spatial maps\", fillcolor=\"#EFEFEF\"]
    Curves [label=\"Time series + CRI\", fillcolor=\"#EFEFEF\"]
    Paper [label=\"Manuscript + PPT\", fillcolor=\"#EFEFEF\"]
  }
  eBird->Dedup; BW->Dedup; iNat->Dedup
  Dedup->Grid; Dedup->Blocks; Dedup->Cov
  Grid->Fit; Blocks->Fit; Cov->Fit
  Fit->Diag; Fit->Psi
  Psi->Tax; Psi->Phy; Psi->Fun; Psi->Beta
  Tax->VIF; Phy->VIF; Fun->VIF; Beta->VIF
  VIF->Brms->Imp
  Tax->Maps; Phy->Maps; Fun->Maps; Beta->Maps
  Tax->Curves; Phy->Curves; Fun->Curves
  Maps->Paper; Curves->Paper; Imp->Paper
}"
svg <- DiagrammeRsvg::export_svg(DiagrammeR::grViz(dot))
rsvg::rsvg_pdf(charToRaw(svg), file=v3_file("figures","fig_v3_workflow_pipeline","pdf"), width=1800)
rsvg::rsvg_png(charToRaw(svg), file=v3_file("figures","fig_v3_workflow_pipeline","png"), width=2400)
log_time("20","done")
