#!/usr/bin/env Rscript
## 00_config.R  —  v3 集中配置
##
## 所有脚本 source 此文件获取路径和参数。不执行任何计算。

# ── 项目根目录（自动检测）─────────────────────────────────────────────
PROJECT_ROOT <- Sys.getenv(
  "BIRD_PROJECT_ROOT",
  if (dir.exists(file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis"))) {
    file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis")
  } else {
    getwd()
  }
)
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, mustWork = FALSE)

# ── 目录布局 ──────────────────────────────────────────────────────────
DIRS <- list(
  code         = file.path(PROJECT_ROOT, "code_v3"),
  data_raw     = file.path(PROJECT_ROOT, "data"),
  derived      = file.path(PROJECT_ROOT, "data", "derived_v3"),
  results      = file.path(PROJECT_ROOT, "results_v3"),
  figures      = file.path(PROJECT_ROOT, "figures_v3"),
  logs         = file.path(PROJECT_ROOT, "logs_v3"),
  v2_derived   = file.path(PROJECT_ROOT, "data", "derived_v2"),
  v2_results   = file.path(PROJECT_ROOT, "results_v2")
)

# ── 时间范围 ──────────────────────────────────────────────────────────
ANALYSIS_YR_LO   <- 2000L
ANALYSIS_YR_HI   <- 2024L
PERIOD_LENGTH    <- 5L          # 5 年一个 primary period
BREEDING_MONTHS  <- 4:8         # 4-8 月（可配置）

# ── 空间网格 ──────────────────────────────────────────────────────────
GRID_SIZE_KM     <- as.integer(Sys.getenv("V3_GRID_SIZE_KM", "100")); if (length(GRID_SIZE_KM) == 0) GRID_SIZE_KM <- 100L
MIN_VISITS       <- as.integer(Sys.getenv("V3_MIN_VISITS", "3"));     if (length(MIN_VISITS) == 0) MIN_VISITS <- 3L
MIN_SPECIES_GRID <- as.integer(Sys.getenv("V3_MIN_SPECIES_GRID", "5")); if (length(MIN_SPECIES_GRID) == 0) MIN_SPECIES_GRID <- 5L

# ── 运行标签 ──────────────────────────────────────────────────────────
GRID_TAG <- if (GRID_SIZE_KM != 100L) paste0("_", GRID_SIZE_KM, "km") else ""
RUN_LABEL <- Sys.getenv("V3_RUN_LABEL",
  paste0("v3_full_200sp_ar1_spatial", GRID_TAG))
PILOT_LABEL <- Sys.getenv("V3_PILOT_LABEL",
  paste0("v3_pilot_60sp_ar1_spatial", GRID_TAG))

# ── stMsPGOcc 模型参数 ──────────────────────────────────────────────
# Pilot（快速验证）
PILOT_N_BATCH   <- 225L
PILOT_N_BURN    <- 2000L
PILOT_N_THIN    <- 1L
PILOT_N_CHAINS  <- 4L

# Full（正式运行）
FULL_N_BATCH   <- 400L
FULL_N_BURN    <- 5000L
FULL_N_THIN    <- 2L
FULL_N_CHAINS  <- 4L

# 空间参数（stMsPGOcc 新增）
N_NEIGHBORS     <- 5L
COV_MODEL       <- "exponential"   # "exponential", "spherical", "gaussian", "matern"
N_NEIGHBORS_HO  <- 5L              # hold-out 邻居数
N_NEAREST_DIAG  <- 100L            # 种级诊断抽样数

# MCMC 收敛阈值
RHAT_THRESHOLD  <- 1.05
ESS_THRESHOLD   <- 200L

# ── 后处理参数 ──────────────────────────────────────────────────────
PSI_MAX_DRAWS    <- 400L        # psi 后验最大抽取数
POST_DRAWS_USE   <- 300L        # 多样性计算使用 draws 数
TREND_DRAWS      <- 400L        # 趋势计算 draws 数

# ── brms 参数 ──────────────────────────────────────────────────────
BRMS_ITER        <- 4000L
BRMS_WARMUP      <- 2000L
BRMS_CHAINS      <- 4L
BRMS_ADAPT_DELTA <- 0.99
BRMS_MAX_TREED   <- 15L
BRMS_SEED        <- 2024L

# ── 性状参数 ──────────────────────────────────────────────────────
# 基础性状（来自 table_species_traits_imputed.csv / AVONET）
TRAIT_VARS_BASIC <- c(
  "body_mass_g", "clutch_size", "longevity_y",
  "maturity_y", "avonet_hwi", "avonet_range_size"
)

# 扩展性状（v3 新增）
TRAIT_VARS_EXTENDED <- c(
  "diet_specialization",    # Morelli et al. 2021, 基于 EltonTraits
  "habitat_breadth"         # IUCN Red List, 使用栖息地数量
)

# 性状回归使用的 z-score 化变量名（14_species_trait_regression.R 实际公式使用）
TRAIT_VARS_REGRESSION_Z <- c(
  "z_body_mass", "z_hwi", "z_range_size", "z_clutch_size",
  "z_diet_specialization", "z_habitat_breadth"
)

# 所有性状列（原始名，用于 trait_extended_v3.rds 的列选择）
TRAIT_VARS_ALL <- c(TRAIT_VARS_BASIC, TRAIT_VARS_EXTENDED)

# FIX #13: 性状变换类型标注
# log10 变换适用于右偏大范围变量（体质量、分布范围等）
# 不需要 log 变换的变量（已在合理范围内或为指数/计数型）
TRAIT_VARS_LOG10 <- c("body_mass_g", "clutch_size", "longevity_y",
                       "maturity_y", "avonet_range_size")
TRAIT_VARS_NO_LOG <- c("avonet_hwi", "diet_specialization", "habitat_breadth")

# ── 外部数据路径 ─────────────────────────────────────────────────────
# AVONET
AVONET_PATH <- Sys.getenv(
  "V3_AVONET_PATH",
  file.path("~", "lucc", "AVONET1_BirdLife.csv")
)

# EltonTraits 1.0 (Wilman et al. 2014)
ELTONTRAITS_PATH <- Sys.getenv(
  "V3_ELTONTRAITS_PATH",
  file.path("~", "lucc", "EltonTraits", "BirdFuncDat.txt")
)

# v1 性状插补表
TRAIT_IMPUTED_PATH <- Sys.getenv(
  "V3_TRAIT_IMPUTED_PATH",
  file.path(PROJECT_ROOT, "bird_grid_community_analysis", "results",
            "table_species_traits_imputed.csv")
)

# IUCN Red List API token
IUCN_API_KEY <- Sys.getenv("IUCN_API_KEY", "")

# ── 随机森林重要性参数 ──────────────────────────────────────────────
RF_NUM_DRAWS    <- 100L         # 后验抽取数（传播不确定性）
RF_NUM_TREES    <- 1000L
RF_SEED         <- 2025L
RF_MIN_NODE     <- 5L
RF_MTRY_FACTOR  <- NULL        # NULL = ranger 默认 sqrt(p)

# ── 驱动变量分组（varpart & RF 共用）─────────────────────────────────
# v3 改进：使用时间变化量作为主驱动因子
# 逻辑：响应变量是群落时间趋势 → 预测变量也应是环境的时间变化

# 后置驱动分析使用的变化量驱动组（09_extended_analyses.R）
DRIVER_GROUPS_TREND <- list(
  climate_change   = c("delta_t_mean", "delta_t_extreme"),     # CRU TS 温度变化
  landuse_change   = c("delta_forest", "delta_impervious",     # CLCD 土地覆盖变化
                        "delta_natural"),
  human_change     = c("delta_hfi"),                           # HFI 变化
  baseline         = c("bio11", "elev_mean",                   # 静态基线控制
                        "centroid_lon", "centroid_lat")
)

# 敏感性分析：使用标准化气候变化量
DRIVER_GROUPS_TREND_SENSITIVITY <- list(
  climate_change   = c("delta_t_std"),                         # 标准化温度变化
  landuse_change   = c("delta_forest", "delta_impervious",
                        "delta_natural"),
  human_change     = c("delta_hfi"),
  baseline         = c("bio11", "elev_mean",
                        "centroid_lon", "centroid_lat")
)

# stMsPGOcc 占有率模型使用的静态协变量组（保持不变）
DRIVER_GROUPS <- list(
  climate      = c("bio4", "bio7", "bio11", "bio13"),
  topo_habitat = c("elev_mean", "elev_sd", "texture_shannon",
                    "habitat_diversity_shannon"),
  human        = c("hfi_mean", "landcover_built", "landcover_cropland"),
  space        = c("centroid_lon", "centroid_lat")
)

# 变化量驱动变量的标签映射（用于图表）
DRIVER_TREND_LABELS <- c(
  delta_t_mean     = "Mean temperature change",
  delta_t_extreme  = "Extreme temperature change",
  delta_t_std      = "Standardized temperature change",
  delta_forest     = "Forest cover change",
  delta_impervious = "Impervious surface change",
  delta_natural    = "Natural land cover change",
  delta_hfi        = "Human footprint change",
  bio11            = "Mean temperature (baseline)",
  elev_mean        = "Elevation (baseline)",
  centroid_lon     = "Longitude",
  centroid_lat     = "Latitude"
)

# 变化量驱动组的标签映射
DRIVER_GROUP_TREND_LABELS <- c(
  climate_change = "Climate change",
  landuse_change = "Land use change",
  human_change   = "Human pressure change",
  baseline       = "Spatial baseline"
)

# ── 外部数据路径（新增）─────────────────────────────────────────────
# CRU TS
CRU_DIR <- Sys.getenv("V3_CRU_DIR",
  file.path(PROJECT_ROOT, "data", "external", "cru_ts"))

# CLCD
CLCD_DIR <- Sys.getenv("V3_CLCD_DIR",
  file.path(PROJECT_ROOT, "data", "external", "clcd"))

# ── Nature/Science 图表规范 ──────────────────────────────────────────
NATURE_FONT      <- "Arial"
NATURE_PT        <- 7.5
NATURE_WIDTH_S   <- 89    # mm, 单栏
NATURE_WIDTH_M   <- 120   # mm, 1.5 栏
NATURE_WIDTH_L   <- 183   # mm, 双栏
NATURE_ACCENT    <- "#0E5A78"
NATURE_DPI       <- 300L
NATURE_LINE_AXIS <- 0.25
NATURE_LINE_TICK <- 0.15
NATURE_LINE_MAP  <- 0.12

# ── 检测协变量（v3 修正 + FIX #4）─────────────────────────────────────
# 移除 log_observers（与 log_events 共线）
# FIX #4: has_duration 信息量极低（仅反映数据来源差异）
# 改善方案：log_duration（有值时）+ has_duration（缺失标志变量）
# has_duration 作为 missingness indicator：有 duration → 1，缺失 → 0
# 缺失 duration 时 log_duration 用 0 填补（has_duration=0 表示该值无效）
DET_FORMULA_STR <- "~ log_events + log_duration + has_duration"

# ── 占有率协变量 ──────────────────────────────────────────────────────
OCC_FORMULA_STR <- "~ bio4 + bio7 + bio11 + bio13 + elev_mean + elev_sd +
  texture_shannon + habitat_diversity_shannon + hfi_mean +
  landcover_built + landcover_cropland + centroid_lon + centroid_lat +
  year_scaled"

# ── 确保 v3 输出目录存在 ──────────────────────────────────────────────
ensure_v3_dirs <- function() {
  for (d in DIRS[c("derived", "results", "figures", "logs")]) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  DIRS
}

# ── v3 文件路径辅助由 utils_paths.R 提供 ──────────────────────────────
# v3_file(), v2_file(), project_file(), log_path() 均定义在 utils_paths.R
# 00_config.R 只提供 DIRS 布局，具体路径函数统一由 utils_paths.R 管理

message("[00_config] loaded  —  project: ", PROJECT_ROOT)
