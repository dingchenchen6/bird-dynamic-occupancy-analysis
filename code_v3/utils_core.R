#!/usr/bin/env Rscript
## utils_core.R  —  v3 通用函数库
##
## 提供后验汇总、标准化、数组填充、检查点等通用工具

# ── 加载依赖 ──────────────────────────────────────────────────────────
# 使用稳健路径检测，避免 sys.frame() 在 Rscript/RStudio 下的脆弱性
{
  .root <- Sys.getenv(
    "BIRD_PROJECT_ROOT",
    if (dir.exists(file.path("~", "Documents", "New project",
                             "bird_dynamic_occupancy_analysis"))) {
      file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis")
    } else getwd()
  )
  source(file.path(.root, "code_v3", "utils_paths.R"))
  rm(.root)
}

# ── 后验汇总（统一版本，v2 在 04/04b 重复定义）───────────────────────
#' @param x 后验样本矩阵 (draws × param) 或数值向量
#' @param probs 分位数概率
#' @return tibble: mean, sd, q025, median, q975, rhat (如可用)
summarise_post <- function(x, probs = c(0.025, 0.5, 0.975)) {
  if (is.array(x) && length(dim(x)) == 3) {
    d <- dim(x)
    x <- matrix(x, nrow = d[1] * d[3], ncol = d[2])
  }
  if (is.matrix(x) || is.array(x)) {
    # 矩阵：逐列汇总
    qs <- apply(x, 2, quantile, probs = probs, na.rm = TRUE)
    tibble(
      mean   = colMeans(x, na.rm = TRUE),
      sd     = matrixStats::colSds(x, na.rm = TRUE),
      q025   = qs["2.5%", ],
      median = qs["50%", ],
      q975   = qs["97.5%", ]
    )
  } else {
    # 向量
    qs <- quantile(x, probs = probs, na.rm = TRUE)
    tibble(
      mean   = mean(x, na.rm = TRUE),
      sd     = sd(x, na.rm = TRUE),
      q025   = qs["2.5%"],
      median = qs["50%"],
      q975   = qs["97.5%"]
    )
  }
}

# ── Z-score 标准化 ───────────────────────────────────────────────────
standardize <- function(x, center = TRUE, scale = TRUE) {
  if (isTRUE(center)) x <- x - mean(x, na.rm = TRUE)
  if (isTRUE(scale))  x <- x / sd(x, na.rm = TRUE)
  x
}

# ── 填充数组中的 NA（v2 中反复重写）───────────────────────────────────
#' 用列中位数填充矩阵/数组的 NA
fill_arr_na <- function(mat, margin = 2L) {
  if (is.matrix(mat)) {
    if (margin == 2L) {
      for (j in seq_len(ncol(mat))) {
        miss <- is.na(mat[, j])
        if (any(miss)) mat[miss, j] <- median(mat[, j], na.rm = TRUE)
      }
    } else {
      for (i in seq_len(nrow(mat))) {
        miss <- is.na(mat[i, ])
        if (any(miss)) mat[i, miss] <- median(mat[i, ], na.rm = TRUE)
      }
    }
  } else if (is.array(mat) && length(dim(mat)) == 3L) {
    # 3D 数组: 逐 3rd-dim 切片填充
    for (k in seq_len(dim(mat)[3L])) {
      mat[, , k] <- fill_arr_na(mat[, , k], margin = margin)
    }
  }
  mat
}

# ── 检查点保存 ────────────────────────────────────────────────────────
#' 安全保存 RDS 检查点
checkpoint_save <- function(obj, name, subdir = "derived") {
  path <- v3_file(subdir, name, "rds")
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(obj, path)
  message(sprintf("[checkpoint] saved %s → %s", name, path))
  invisible(path)
}

# ── 安全读取 CSV ──────────────────────────────────────────────────────
read_csv_safe <- function(path, ...) {
  if (!file.exists(path)) {
    warning("File not found: ", path, call. = FALSE)
    return(NULL)
  }
  readr::read_csv(path, show_col_types = FALSE, ...)
}

# ── 安全读取 RDS ──────────────────────────────────────────────────────
safe_read <- function(path, verbose = TRUE) {
  if (!file.exists(path)) {
    if (verbose) warning("File not found: ", path, call. = FALSE)
    return(NULL)
  }
  readRDS(path)
}

# ── 时间戳日志 ────────────────────────────────────────────────────────
log_time <- function(stage, msg = "") {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s %s", ts, stage, msg))
}

# ── 列中位数填充 tibble ──────────────────────────────────────────────
fill_median_tibble <- function(tbl, cols) {
  for (cc in cols) {
    if (cc %in% names(tbl)) {
      v <- as.numeric(tbl[[cc]])
      if (any(is.na(v))) v[is.na(v)] <- median(v, na.rm = TRUE)
      tbl[[cc]] <- v
    }
  }
  tbl
}

message("[utils_core] loaded")
