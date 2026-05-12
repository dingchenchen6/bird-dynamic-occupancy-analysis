#!/usr/bin/env Rscript
## utils_paths.R  —  v3 路径工具模块
##
## 集中管理 %||% 操作符和路径工具函数

# ── 统一 %||% 操作符（v2 中在3个文件重复定义）────────────────────────
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ── 项目根目录自动检测 ────────────────────────────────────────────────
detect_project_root <- function() {
  # 1. 环境变量
  env_val <- Sys.getenv("BIRD_PROJECT_ROOT", "")
  if (nzchar(env_val) && dir.exists(env_val)) return(normalizePath(env_val))

  # 2. 从当前工作目录向上查找
  cwd <- getwd()
  for (i in 0:3) {
    candidate <- if (i == 0) cwd else dirname(candidate)
    marker <- file.path(candidate, "code_v3")
    if (dir.exists(marker)) return(normalizePath(candidate))
  }

  # 3. 默认路径
  default <- "~/Documents/New project/bird_dynamic_occupancy_analysis"
  if (dir.exists(default)) return(normalizePath(default))

  stop("Cannot detect project root. Set BIRD_PROJECT_ROOT environment variable.")
}

# ── 加载配置（若尚未加载）────────────────────────────────────────────
load_config_if_needed <- function() {
  if (!exists("DIRS", envir = .GlobalEnv, inherits = FALSE)) {
    source(file.path(detect_project_root(), "code_v3", "00_config.R"))
  }
}

# ── v3 文件路径 ──────────────────────────────────────────────────────
v3_file <- function(subdir, stem, ext = "csv") {
  load_config_if_needed()
  base <- DIRS[[subdir]]
  if (is.null(base)) stop("Unknown subdir: ", subdir, call. = FALSE)
  file.path(base, paste0(stem, ".", ext))
}

# ── 兼容 v2 路径（读取旧数据时用）─────────────────────────────────────
v2_file <- function(subdir, stem, ext = "csv") {
  load_config_if_needed()
  v2_map <- list(
    derived  = DIRS$v2_derived,
    results  = DIRS$v2_results,
    figures  = file.path(PROJECT_ROOT, "figures_v2"),
    logs     = file.path(PROJECT_ROOT, "logs_v2")
  )
  base <- v2_map[[subdir]]
  if (is.null(base)) stop("Unknown v2 subdir: ", subdir, call. = FALSE)
  file.path(base, paste0(stem, ".", ext))
}

# ── 项目数据路径 ──────────────────────────────────────────────────────
project_file <- function(...) {
  load_config_if_needed()
  file.path(PROJECT_ROOT, ...)
}

# ── 日志输出路径 ──────────────────────────────────────────────────────
log_path <- function(script_name) {
  load_config_if_needed()
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  v3_file("logs", paste0(script_name, "_", ts), "log")
}

# ── 统一 source 辅助（其他 utils 模块调用）─────────────────────────────
#' 安全加载 code_v3/ 下的 R 模块
#' @param module 文件名（不含路径），如 "utils_core.R"
source_v3_module <- function(module) {
  root <- detect_project_root()
  path <- file.path(root, "code_v3", module)
  if (!file.exists(path)) stop("Module not found: ", path, call. = FALSE)
  source(path)
}

# ── 带网格标签回退的文件加载 ──────────────────────────────────────────
#' 优先加载带 GRID_TAG 的文件，回退到无标签版本
#' @param subdir 子目录名（如 "derived"）
#' @param stem 文件名主干（如 "survey_history_v3"）
#' @param ext 扩展名
v3_file_tagged <- function(subdir, stem, ext = "rds") {
  # 先尝试带标签版本（在 stem 中插入 GRID_TAG）
  if (exists("GRID_TAG", envir = .GlobalEnv) && nzchar(GRID_TAG)) {
    tagged_stem <- sub("_v3$", paste0(GRID_TAG, "_v3"), stem)
    tagged_path <- v3_file(subdir, tagged_stem, ext)
    if (file.exists(tagged_path)) return(tagged_path)
  }
  # 回退到原始路径
  v3_file(subdir, stem, ext)
}

message("[utils_paths] loaded")
