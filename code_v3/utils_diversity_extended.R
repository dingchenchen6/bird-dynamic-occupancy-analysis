#!/usr/bin/env Rscript
## utils_diversity_extended.R — 扩展多样性计算工具
##
## 科学问题 / Scientific question:
## 在 v3 基础上引入 fundiversity 包计算更多功能多样性指数，
## 并支持 McTavish CLOOTL 系统发育树。
##
## 主要扩展 / Key extensions:
##   1. fundiversity: FRic, FDiv, FEve, FDis, Rao's Q
##   2. Trait PCA 降维（FRic 需要 n_axes < n_species）
##   3. McTavish 树支持（pd_prob_mctavish, mpd_prob_mctavish）
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# ── 功能多样性扩展（fundiversity）───────────────────────────────────────
#' 使用 fundiversity 计算功能多样性指数（概率加权）
#'
#' @param psi 数值向量：物种占有率后验（概率）
#' @param trait_mat 标准化后的性状矩阵（species × traits）
#' @param n_pca_axes PCA 降维轴数（FRic 需要 n_axes < n_species）
#' @param eps 浮点精度
#' @return list: fric_prob, fdiv_fund, feve_fund, fdis_prob, fmpd_prob, raoq_fund,
#'               cwm_pc1, cwm_pc2
#'
#' 说明 / Note:
#' fundiversity 需要 site × species 矩阵输入。
#' 单 site 场景下将 psi 包装为 1 行矩阵。
#' FRic 计算前对性状矩阵做 PCA，避免维度灾难。
div_functional_fundiversity <- function(psi, trait_mat,
                                         n_pca_axes = 5L,
                                         eps = 1e-12) {
  ok <- psi > eps
  S <- sum(ok)

  if (S < 3 || is.null(trait_mat) || nrow(trait_mat) != length(psi)) {
    return(list(
      fric_prob = NA_real_, fdiv_fund = NA_real_, feve_fund = NA_real_,
      fdis_prob = NA_real_, fmpd_prob = NA_real_, raoq_fund = NA_real_,
      cwm_pc1 = NA_real_, cwm_pc2 = NA_real_
    ))
  }

  psi_ok   <- psi[ok]
  trait_ok <- trait_mat[ok, , drop = FALSE]
  w <- psi_ok / sum(psi_ok)

  # PCA 降维（FRic 要求 n_axes < n_species）
  trait_pca <- prcomp(trait_ok, scale. = FALSE, center = TRUE)
  n_axes_use <- min(n_pca_axes, S - 1, ncol(trait_ok))
  trait_axes <- as.data.frame(trait_pca$x[, seq_len(n_axes_use), drop = FALSE])
  rownames(trait_axes) <- rownames(trait_ok)

  # CWM on first 2 PCA axes
  cwm_pc1 <- sum(w * trait_axes[, 1])
  cwm_pc2 <- if (ncol(trait_axes) >= 2) sum(w * trait_axes[, 2]) else NA_real_

  # fundiversity 输入：site × species 矩阵（单 site = 1 行）
  # FIX: 使用二值 presence/absence（psi > 0.5），避免所有物种被视为存在
  sp_com <- matrix(as.numeric(psi_ok > 0.5), nrow = 1,
                    dimnames = list("site1", rownames(trait_axes)))

  # 计算各指数（tryCatch 保护，防止 convex hull 失败）
  # FIX: 移除 stand = TRUE，避免除以全局凸包导致恒为 1
  fric_res <- tryCatch(
    fundiversity::fd_fric(traits = trait_axes, sp_com = sp_com, stand = FALSE),
    error = function(e) { warning("[fundiv] fd_fric failed: ", e$message); data.frame(site = "site1", FRic = NA_real_) }
  )
  fdiv_res <- tryCatch(
    fundiversity::fd_fdiv(traits = trait_axes, sp_com = sp_com),
    error = function(e) { warning("[fundiv] fd_fdiv failed: ", e$message); data.frame(site = "site1", FDiv = NA_real_) }
  )
  feve_res <- tryCatch(
    fundiversity::fd_feve(traits = trait_axes, sp_com = sp_com),
    error = function(e) { warning("[fundiv] fd_feve failed: ", e$message); data.frame(site = "site1", FEve = NA_real_) }
  )
  fdis_res <- tryCatch(
    fundiversity::fd_fdis(traits = trait_axes, sp_com = sp_com),
    error = function(e) { warning("[fundiv] fd_fdis failed: ", e$message); data.frame(site = "site1", FDis = NA_real_) }
  )
  raoq_res <- tryCatch(
    fundiversity::fd_raoq(traits = trait_axes, sp_com = sp_com),
    error = function(e) { warning("[fundiv] fd_raoq failed: ", e$message); data.frame(site = "site1", RaoQ = NA_real_) }
  )

  # fmpd: functional mean pairwise distance（trait-space MPD 类比）
  fmpd_prob <- NA_real_
  if (nrow(trait_axes) >= 2) {
    dist_mat <- as.matrix(dist(trait_axes))
    w_mat <- outer(w, w)
    diag(w_mat) <- 0
    if (sum(w_mat) > eps) {
      fmpd_prob <- sum(dist_mat * w_mat) / sum(w_mat)
    }
  }

  list(
    fric_prob = as.numeric(fric_res$FRic[1]),
    fdiv_fund = as.numeric(fdiv_res$FDiv[1]),
    feve_fund = as.numeric(feve_res$FEve[1]),
    fdis_prob = as.numeric(fdis_res$FDis[1]),
    fmpd_prob = fmpd_prob,
    raoq_fund = as.numeric(raoq_res$Q[1]),
    cwm_pc1   = cwm_pc1,
    cwm_pc2   = cwm_pc2
  )
}

# ── 系统发育多样性（McTavish CLOOTL 支持）───────────────────────────────
#' 概率加权 Faith's PD 和 MPD（支持任意 phylo 对象）
#'
#' 与 div_phylogenetic 相同，但接受外部 phylo 参数，
#' 用于同时计算 Jetz 和 McTavish 两棵树。
div_phylogenetic_generic <- function(psi, phylo, tip_order = NULL, eps = 1e-12) {
  if (is.null(phylo)) {
    return(list(pd_prob = NA_real_, mpd_prob = NA_real_))
  }

  psi <- pmax(psi, 0)
  ok <- psi > eps

  if (sum(ok) < 2) {
    return(list(pd_prob = NA_real_, mpd_prob = NA_real_))
  }

  if (!is.null(tip_order)) {
    names(psi) <- tip_order
  }
  present <- names(psi)[ok]
  present <- intersect(present, phylo$tip.label)
  if (length(present) < 2) {
    return(list(pd_prob = NA_real_, mpd_prob = NA_real_))
  }

  pruned <- ape::drop.tip(phylo, setdiff(phylo$tip.label, present))
  w <- psi[pruned$tip.label] / sum(psi[pruned$tip.label])

  samp <- data.frame(matrix(1, nrow = 1, ncol = length(pruned$tip.label)))
  colnames(samp) <- pruned$tip.label
  rownames(samp) <- "site1"
  pd_raw <- picante::pd(samp, pruned, include.root = FALSE)$PD
  pd_prob <- pd_raw * sum(psi[ok])

  cophen <- ape::cophenetic.phylo(pruned)
  w_mat <- outer(w, w)
  diag(w_mat) <- 0
  mpd_prob <- if (sum(w_mat) > 0) sum(cophen * w_mat) / sum(w_mat) else NA_real_

  list(pd_prob = pd_prob, mpd_prob = mpd_prob)
}

message("[utils_diversity_extended] loaded — fundiversity + McTavish support")
