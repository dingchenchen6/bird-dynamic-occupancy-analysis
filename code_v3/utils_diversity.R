#!/usr/bin/env Rscript
## utils_diversity.R  —  v3 多样性计算工具
##
## 修复 v2 中 log10 bug（if_else(.x > 0, log10(.x), .x) 双峰分布）
## 统一多样性计算逻辑

# ── 加载依赖 ──────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
})
{
  .root <- Sys.getenv("BIRD_PROJECT_ROOT",
    if (dir.exists(file.path("~", "Documents", "New project",
                             "bird_dynamic_occupancy_analysis")))
      file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis")
    else getwd())
  source(file.path(.root, "code_v3", "utils_paths.R"))
  source(file.path(.root, "code_v3", "00_config.R"))
  rm(.root)
}

# ── 功能多样性计算 ──────────────────────────────────────────────────
#' 计算群落加权均值/标准差、性状体积、Rao Q
#'
#' @param psi 数值向量：物种占有率后验（概率）
#' @param trait_mat 标准化后的性状矩阵（species × traits）
#' @param eps 浮点精度
#' @return list: cwm, cwsd, trait_volume, trait_dispersion, rao_q
div_functional <- function(psi, trait_mat, eps = 1e-12) {
  ok <- psi > eps
  n_trait <- ncol(trait_mat)

  # 边界情况
  if (!any(ok) || is.null(trait_mat) || nrow(trait_mat) != length(psi)) {
    return(list(
      cwm = setNames(rep(NA_real_, n_trait), colnames(trait_mat)),
      cwsd = setNames(rep(NA_real_, n_trait), colnames(trait_mat)),
      trait_volume = NA_real_,
      trait_dispersion = NA_real_,
      rao_q = NA_real_
    ))
  }

  psi_ok  <- psi[ok]
  trait_ok <- trait_mat[ok, , drop = FALSE]
  w <- psi_ok / sum(psi_ok)

  # CWM
  cwm <- colSums(w * trait_ok)

  # CWSD
  cwsd <- sqrt(colSums(w * (trait_ok - matrix(cwm, nrow = nrow(trait_ok),
                                                 ncol = ncol(trait_ok),
                                                 byrow = TRUE))^2))

  # Trait volume（概率加权多维 sd）
  cov_mat <- matrix(0, n_trait, n_trait)
  for (i in seq_len(nrow(trait_ok))) {
    diff_i <- trait_ok[i, ] - cwm
    cov_mat <- cov_mat + w[i] * outer(diff_i, diff_i)
  }
  trait_volume <- sqrt(sum(diag(cov_mat)))

  # Trait dispersion（加权重心距离）
  trait_dispersion <- sum(w * sqrt(rowSums((trait_ok - matrix(cwm, nrow = nrow(trait_ok),
                                                                ncol = ncol(trait_ok),
                                                                byrow = TRUE))^2)))

  # Rao's Q（需要距离矩阵）
  dist_mat <- as.matrix(stats::dist(trait_ok))
  rao_q <- sum(outer(w, w) * dist_mat)

  list(cwm = cwm, cwsd = cwsd,
       trait_volume = trait_volume,
       trait_dispersion = trait_dispersion,
       rao_q = rao_q)
}

# ── 分类学多样性 ──────────────────────────────────────────────────────
#' 计算占有率校正的丰富度、Shannon、逆 Simpson
#'
#' @param psi 数值向量：物种占有率
#' @return list: richness, shannon, inv_simpson
div_taxonomic <- function(psi, eps = 1e-12) {
  psi <- pmax(psi, 0)
  # FIX #12: 占域校正丰富度 = sum(psi)，即期望物种数
  # 旧版 sum(psi > eps) 几乎等于物种池大小（因为 eps=1e-12 极小）
  richness <- sum(psi)
  p <- psi / sum(psi)

  shannon <- -sum(p[p > eps] * log(p[p > eps]))
  inv_simpson <- 1 / sum(p^2)

  list(richness = richness, shannon = shannon, inv_simpson = inv_simpson)
}

# ── 系统发育多样性 ────────────────────────────────────────────────────
#' 概率加权 Faith's PD 和 MPD
#'
#' @param psi 数值向量：物种占有率
#' @param phylo 系统发育树（ape::phylo）
#' @param tip_order 物种顺序（与 psi 对齐）
#' @return list: pd_prob, mpd_prob
div_phylogenetic <- function(psi, phylo, tip_order = NULL, eps = 1e-12) {
  if (is.null(phylo)) {
    return(list(pd_prob = NA_real_, mpd_prob = NA_real_))
  }

  psi <- pmax(psi, 0)
  ok <- psi > eps

  if (sum(ok) < 2) {
    return(list(pd_prob = NA_real_, mpd_prob = NA_real_))
  }

  # 对齐 tip 与 psi 顺序
  if (!is.null(tip_order)) {
    names(psi) <- tip_order
  }
  present <- names(psi)[ok]

  # 检查 tip 是否在树中
  present <- intersect(present, phylo$tip.label)
  if (length(present) < 2) {
    return(list(pd_prob = NA_real_, mpd_prob = NA_real_))
  }

  # PD: 加权分支长度
  pruned <- ape::drop.tip(phylo, setdiff(phylo$tip.label, present))
  w <- psi[pruned$tip.label] / sum(psi[pruned$tip.label])
  # picante::pd 需要样地 × 物种的 0/1 矩阵（行=样地，列=物种名）
  samp <- data.frame(matrix(1, nrow = 1, ncol = length(pruned$tip.label)))
  colnames(samp) <- pruned$tip.label
  rownames(samp) <- "site1"
  pd_raw <- picante::pd(samp, pruned, include.root = FALSE)$PD
  pd_prob <- pd_raw * sum(psi[ok])  # 按总占有率缩放

  # MPD: 加权平均系统发育距离
  cophen <- cophenetic.phylo(pruned)
  w_mat <- outer(w, w)
  diag(w_mat) <- 0
  mpd_prob <- sum(cophen * w_mat) / sum(w_mat)

  list(pd_prob = pd_prob, mpd_prob = mpd_prob)
}

# ── 群落时间动态指标（codyn 风格）─────────────────────────────────────
## 来源：codyn R 包 (Hallett et al. 2020, Methods Ecol Evol)
## 核心思路：用 ψ 后验占有率替代传统多度，作为群落组成向量
## 三个关键指标：同步性、方差比、时间周转率

#' Loreau 同步性指数
#'
#' 度量物种动态的同步程度。值域 [0,1]：
#'   1 = 完全同步（所有物种同步涨落）
#'   0 = 完全异步（物种涨落相互抵消）
#' 参考：Loreau & de Mazancourt 2008, Am Nat
#'
#' @param comm_matrix 群落矩阵 (n_species × n_periods)，值为占有率或相对多度
#' @param eps 浮点精度
#' @return numeric: 同步性指数 [0,1]
div_synchrony <- function(comm_matrix, eps = 1e-12) {
  if (is.null(comm_matrix) || nrow(comm_matrix) < 2 || ncol(comm_matrix) < 2) return(NA_real_)

  # 去除方差为零的物种（恒定占有率，无时间动态）
  sp_var <- apply(comm_matrix, 1, var, na.rm = TRUE)
  keep <- sp_var > eps
  if (sum(keep) < 2) return(NA_real_)
  comm_matrix <- comm_matrix[keep, , drop = FALSE]

  # Loreau de Mazancourt synchrony: φ = σ²_total / (Σ σᵢ)²
  total_ts  <- colSums(comm_matrix)
  var_total <- var(total_ts, na.rm = TRUE)
  sum_var_sp <- sum(apply(comm_matrix, 1, var, na.rm = TRUE))

  if (sum_var_sp < eps) return(NA_real_)
  synchrony <- var_total / (sum_var_sp)^2

  # 裁剪到 [0, 1] 范围（数值误差可能导致微小越界）
  min(1, max(0, synchrony))
}

#' Schluter 方差比 (Variance Ratio)
#'
#' 检验物种间关联性：
#'   VR > 1 → 物种正关联（趋同波动）
#'   VR < 1 → 物种负关联（趋异波动）
#'   VR = 1 → 物种独立波动
#' 参考：Schluter 1984, Ecology
#'
#' @param comm_matrix 群落矩阵 (n_species × n_periods)
#' @param eps 浮点精度
#' @return numeric: 方差比 VR
div_variance_ratio <- function(comm_matrix, eps = 1e-12) {
  if (is.null(comm_matrix) || nrow(comm_matrix) < 2 || ncol(comm_matrix) < 2) return(NA_real_)

  # 去除零方差物种
  sp_var <- apply(comm_matrix, 1, var, na.rm = TRUE)
  keep <- sp_var > eps
  if (sum(keep) < 2) return(NA_real_)
  comm_matrix <- comm_matrix[keep, , drop = FALSE]

  # VR = σ²_total / Σ σᵢ²
  total_ts  <- colSums(comm_matrix)
  var_total <- var(total_ts, na.rm = TRUE)
  sum_var_sp <- sum(apply(comm_matrix, 1, var, na.rm = TRUE))

  if (sum_var_sp < eps) return(NA_real_)
  var_total / sum_var_sp
}

#' 时间周转率（codyn 风格平均周转）
#'
#' 计算相邻时间段间的物种增/减/总周转率，取所有相邻时段的均值。
#' 与 Baselga 分解互补：Baselga 分解 Sørensen 距离为周转+嵌套成分，
#' 本指标直接计算物种增减比例。
#'
#' @param comm_matrix 群落矩阵 (n_species × n_periods)
#' @param threshold 占有率阈值，低于此视为不存在（默认 0.1）
#' @return list: turnover_total（总周转）, turnover_gain（增）, turnover_loss（减）
div_temporal_turnover <- function(comm_matrix, threshold = 0.1) {
  if (is.null(comm_matrix) || nrow(comm_matrix) < 2 || ncol(comm_matrix) < 2)
    return(list(turnover_total = NA_real_, turnover_gain = NA_real_, turnover_loss = NA_real_))

  n_periods <- ncol(comm_matrix)
  total_gain <- 0; total_loss <- 0; n_pairs <- 0

  for (t in seq_len(n_periods - 1)) {
    present_t   <- comm_matrix[, t]   > threshold
    present_t1  <- comm_matrix[, t + 1] > threshold

    gain <- sum(!present_t & present_t1)   # t 不在但 t+1 在
    loss <- sum(present_t & !present_t1)   # t 在但 t+1 不在
    pool <- sum(present_t | present_t1)    # 总物种池

    if (pool > 0) {
      total_gain  <- total_gain + gain / pool
      total_loss  <- total_loss + loss / pool
      n_pairs     <- n_pairs + 1
    }
  }

  if (n_pairs == 0)
    return(list(turnover_total = NA_real_, turnover_gain = NA_real_, turnover_loss = NA_real_))

  avg_gain <- total_gain / n_pairs
  avg_loss <- total_loss / n_pairs

  list(
    turnover_total = avg_gain + avg_loss,
    turnover_gain  = avg_gain,
    turnover_loss  = avg_loss
  )
}

# ── 性状矩阵预处理（修复 v2 log10 bug + #13 性状分类变换）────────────
#' FIX #13: 区分需要/不需要 log10 变换的性状
#' - 需 log10：右偏大范围变量（body_mass, range_size 等）
#' - 不需 log10：已在合理范围内的变量（hwi, diet_specialization, habitat_breadth）
#' v2 bug: if_else(.x > 0, log10(.x), .x) 产生双峰分布
#' v3 fix: if_else(.x > 0, log10(.x), NA_real_) 产生正确缺失值
#' v3 fix2: 只对 TRAIT_VARS_LOG10 中的变量做 log10 变换
prepare_trait_matrix <- function(trait_aligned, trait_vars,
                                  log_vars = TRAIT_VARS_LOG10) {
  no_log_vars <- setdiff(trait_vars, log_vars)

  # 需 log 变换的列
  if (length(log_vars) > 0) {
    log_cols <- intersect(log_vars, trait_vars)
    log_part <- trait_aligned |>
      select(any_of(log_cols)) |>
      mutate(across(where(is.numeric),
                    ~ if_else(.x > 0, log10(.x), NA_real_)))
  } else {
    log_part <- tibble(.rows = nrow(trait_aligned))
  }

  # 不需 log 变换的列
  if (length(no_log_vars) > 0) {
    no_log_cols <- intersect(no_log_vars, trait_vars)
    no_log_part <- trait_aligned |>
      select(any_of(no_log_cols))
  } else {
    no_log_part <- tibble(.rows = nrow(trait_aligned))
  }

  # 合并
  trait_mat <- bind_cols(log_part, no_log_part) |>
    select(all_of(trait_vars)) |>
    as.matrix()
  rownames(trait_mat) <- trait_aligned$species

  # NA 用列中位数填补
  for (j in seq_len(ncol(trait_mat))) {
    miss <- is.na(trait_mat[, j])
    if (any(miss)) trait_mat[miss, j] <- median(trait_mat[, j], na.rm = TRUE)
  }

  # 标准化
  trait_mat <- scale(trait_mat)
  attr(trait_mat, "scaled:center") <- NULL
  attr(trait_mat, "scaled:scale")  <- NULL

  trait_mat
}

# ── 功能均匀度 FEve（Villéger et al. 2008, Ecology）───────────────────
#' 概率加权功能均匀度 FEve
#'
#' 度量物种多度在功能空间中最小生成树边上的均匀程度。
#' 值域 [0, 1]：
#'   1 = 多度在功能空间中完全均匀分布
#'   0 = 多度极端集中在少数功能分支
#'
#' 算法（Villéger et al. 2008, Box 1）：
#'   1. 对存在物种计算功能距离矩阵，构建最小生成树（MST）
#'   2. 对 MST 每条边 (i,j)：EW_k = (w_i + w_j) * d_ij / 2
#'   3. 归一化：PEW_k = EW_k / Σ(EW)
#'   4. FEve = 1 - Σ|PEW_k - 1/(S-1)| / (2 * (1 - 1/(S-1)))
#'
#' @param psi 数值向量：物种占有率后验（概率）
#' @param trait_mat 标准化后的性状矩阵（species × traits）
#' @param eps 浮点精度
#' @return numeric: FEve ∈ [0, 1]
div_feve <- function(psi, trait_mat, eps = 1e-12) {
  ok <- psi > eps
  S <- sum(ok)
  if (S < 3 || is.null(trait_mat) || nrow(trait_mat) != length(psi)) return(NA_real_)

  psi_ok   <- psi[ok]
  trait_ok <- trait_mat[ok, , drop = FALSE]
  w <- psi_ok / sum(psi_ok)  # 相对多度

  # 功能距离矩阵
  D <- as.matrix(stats::dist(trait_ok))

  # Prim 最小生成树（不依赖 vegan/igraph）
  n <- S
  in_tree <- rep(FALSE, n)
  min_dist <- rep(Inf, n)
  min_from <- rep(NA_integer_, n)
  min_dist[1] <- 0

  edge_i <- integer(n - 1)
  edge_j <- integer(n - 1)

  for (k in seq_len(n)) {
    # 选择不在树中距离最小的节点
    u <- which.min(ifelse(in_tree, Inf, min_dist))
    if (is.na(u) || min_dist[u] == Inf) break
    in_tree[u] <- TRUE

    if (k > 1) {
      edge_i[k - 1] <- min_from[u]
      edge_j[k - 1] <- u
    }

    # 更新邻居距离
    for (v in seq_len(n)) {
      if (!in_tree[v] && D[u, v] < min_dist[v]) {
        min_dist[v] <- D[u, v]
        min_from[v] <- u
      }
    }
  }

  n_edges <- sum(!is.na(edge_i))
  if (n_edges < 2) return(NA_real_)

  # 计算每条边的 EW 和 PEW
  EW <- numeric(n_edges)
  for (k in seq_len(n_edges)) {
    i <- edge_i[k]; j <- edge_j[k]
    EW[k] <- (w[i] + w[j]) * D[i, j] / 2
  }
  sum_EW <- sum(EW)
  if (sum_EW < eps) return(NA_real_)

  PEW <- EW / sum_EW

  # FEve 公式
  PEW_ref <- 1 / (S - 1)
  FEve <- 1 - sum(abs(PEW - PEW_ref)) / (2 * (1 - PEW_ref))

  max(0, min(1, FEve))
}

# ── 功能分散度 FDiv（Mason et al. 2003; Villéger et al. 2008）────────
#' 概率加权功能分散度 FDiv
#'
#' 度量物种多度在功能空间中偏离重心的程度。
#' 值域 (0, 1]：
#'   高值 = 多度集中在功能空间边缘（功能分散）
#'   低值 = 多度集中在功能空间中心（功能趋同）
#'
#' 算法（Villéger et al. 2008, Box 1）：
#'   1. G = 加权重心
#'   2. dG_i = 物种 i 到 G 的欧氏距离
#'   3. dbar = 加权平均距离
#'   4. dmax = 最大距离
#'   5. FDiv = (sum(w * |dG - dbar|) + dbar) / (sum(w * dG) + dbar)
#'
#' @param psi 数值向量：物种占有率后验（概率）
#' @param trait_mat 标准化后的性状矩阵（species × traits）
#' @param eps 浮点精度
#' @return numeric: FDiv ∈ (0, 1]
div_fdiv <- function(psi, trait_mat, eps = 1e-12) {
  ok <- psi > eps
  if (sum(ok) < 3 || is.null(trait_mat) || nrow(trait_mat) != length(psi)) return(NA_real_)

  psi_ok   <- psi[ok]
  trait_ok <- trait_mat[ok, , drop = FALSE]
  w <- psi_ok / sum(psi_ok)

  # Step 1: 加权重心 G
  G <- colSums(w * trait_ok)

  # Step 2: 各物种到 G 的欧氏距离
  dG <- sqrt(rowSums((trait_ok - matrix(G, nrow = nrow(trait_ok),
                                         ncol = ncol(trait_ok),
                                         byrow = TRUE))^2))

  # Step 3: 加权平均距离
  dbar <- sum(w * dG)

  # Step 4: 最大距离
  dmax <- max(dG)
  if (dmax < eps) return(NA_real_)

  # Step 5: FDiv
  # 分子：加权绝对偏差 + dbar
  # 分母：加权距离和 + dbar
  num <- sum(w * abs(dG - dbar)) + dbar
  den <- sum(w * dG) + dbar

  if (den < eps) return(NA_real_)
  fdiv <- num / den

  max(0, min(1, fdiv))
}

# ── 概率加权 Baselga 分解（Baselga 2010, Glob Ecol Biogeogr）─────────
#' 概率加权时间 β 多样性 Baselga 分解
#'
#' 将 Sørensen 相异分解为周转（Simpson）和嵌套成分：
#'   β_sor = β_sim + β_sne
#'
#' 概率加权版：用 psi 作为"存在概率"替代二值出现/缺失
#'   a = Σ min(psi_1, psi_2)  期望共享物种数
#'   b = max(0, Σ psi_1 - a)  期望仅在时期1的物种数
#'   c = max(0, Σ psi_2 - a)  期望仅在时期2的物种数
#'
#' @param psi_t1 数值向量：时期1的物种占有率
#' @param psi_t2 数值向量：时期2的物种占有率
#' @return list: beta_sor, beta_sim (turnover), beta_sne (nestedness)
div_baselga <- function(psi_t1, psi_t2) {
  psi_t1 <- pmax(psi_t1, 0)
  psi_t2 <- pmax(psi_t2, 0)

  # 概率加权 a, b, c
  a <- sum(pmin(psi_t1, psi_t2))       # 期望共享
  b <- max(0, sum(psi_t1) - a)         # 期望仅在 t1
  c <- max(0, sum(psi_t2) - a)         # 期望仅在 t2

  denom_sor <- 2 * a + b + c
  denom_sim <- a + min(b, c)

  # Sørensen 相异
  beta_sor <- if (denom_sor > 0) (b + c) / denom_sor else NA_real_

  # Simpson 周转
  beta_sim <- if (denom_sim > 0) min(b, c) / denom_sim else NA_real_

  # 嵌套成分
  beta_sne <- if (!is.na(beta_sor) && !is.na(beta_sim)) beta_sor - beta_sim else NA_real_

  # 比例：周转占总 β 的百分比（裁剪到 [0,1]，防止数值误差）
  prop_turnover <- if (!is.na(beta_sor) && beta_sor > 1e-10) {
    min(1, max(0, beta_sim / beta_sor))
  } else NA_real_

  list(
    beta_sor = beta_sor,
    beta_sim = beta_sim,    # 周转成分 (turnover)
    beta_sne = beta_sne,    # 嵌套成分 (nestedness-resultant)
    a = a, b = b, c = c,
    prop_turnover = prop_turnover
  )
}

#' 对所有相邻时期对计算 Baselga 分解
#'
#' @param psi_matrix 群落矩阵 (n_species × n_periods)，值为占有率
#' @return data.frame: period_pair, beta_sor, beta_sim, beta_sne, prop_turnover
div_baselga_temporal <- function(psi_matrix) {
  if (is.null(psi_matrix) || ncol(psi_matrix) < 2)
    return(tibble(period_pair = character(), beta_sor = numeric(),
                  beta_sim = numeric(), beta_sne = numeric(),
                  prop_turnover = numeric()))

  results <- list()
  for (t in seq_len(ncol(psi_matrix) - 1)) {
    bg <- div_baselga(psi_matrix[, t], psi_matrix[, t + 1])
    results[[length(results) + 1]] <- tibble(
      period_pair = paste0("P", t, "_P", t + 1),
      beta_sor = bg$beta_sor,
      beta_sim = bg$beta_sim,
      beta_sne = bg$beta_sne,
      prop_turnover = bg$prop_turnover
    )
  }
  bind_rows(results)
}

# ── 稳健趋势估计（FIX #6: OLS 对端点敏感）──────────────────────────
#' Theil-Sen 斜率估计
#'
#' 所有数据点对的斜率中位数，对异常值和端点不敏感
#' 参考：Theil 1950; Sen 1968, JASA
#'
#' @param y 数值向量（时间序列值）
#' @param x 数值向量（时间索引，默认 1,2,...,n）
#' @return numeric: Theil-Sen 斜率
theil_sen_slope <- function(y, x = seq_along(y)) {
  ok <- !is.na(y) & !is.na(x)
  y <- y[ok]; x <- x[ok]
  n <- length(y)
  if (n < 3) return(NA_real_)

  slopes <- numeric(n * (n - 1) / 2)
  k <- 0
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      k <- k + 1
      slopes[k] <- (y[j] - y[i]) / (x[j] - x[i])
    }
  }
  median(slopes[1:k])
}

#' Mann-Kendall 趋势检验统计量
#'
#' 非参数趋势检验，H0: 无单调趋势
#' S > 0 → 上升趋势；S < 0 → 下降趋势
#'
#' @param y 数值向量
#' @return list: S (统计量), tau (Kendall's tau), p_value (近似 p 值)
mann_kendall_test <- function(y) {
  ok <- !is.na(y)
  y <- y[ok]
  n <- length(y)
  if (n < 4) return(list(S = NA_real_, tau = NA_real_, p_value = NA_real_))

  S <- 0
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      S <- S + sign(y[j] - y[i])
    }
  }

  # Kendall's tau
  tau <- S / (n * (n - 1) / 2)

  # 方差修正（考虑并列值）
  tie_groups <- table(y)
  tie_correction <- sum(tie_groups * (tie_groups - 1) * (2 * tie_groups + 5))
  var_S <- (n * (n - 1) * (2 * n + 5) - tie_correction) / 18

  # Z 近似
  if (S > 0) {
    Z <- (S - 1) / sqrt(max(var_S, 1))
  } else if (S < 0) {
    Z <- (S + 1) / sqrt(max(var_S, 1))
  } else {
    Z <- 0
  }
  p_value <- 2 * pnorm(-abs(Z))

  list(S = S, tau = tau, p_value = p_value)
}

message("[utils_diversity] loaded")
