#!/usr/bin/env Rscript
## utils_spatial.R  —  v3 空间工具模块
##
## 提供测地线距离矩阵计算（修复 v2 用欧氏距离于经纬度的 bug）

# ── 加载依赖 ──────────────────────────────────────────────────────────
{
  .root <- Sys.getenv("BIRD_PROJECT_ROOT",
    if (dir.exists(file.path("~", "Documents", "New project",
                             "bird_dynamic_occupancy_analysis")))
      file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis")
    else getwd())
  source(file.path(.root, "code_v3", "utils_paths.R"))
  rm(.root)
}

# ── 测地线距离矩阵 ──────────────────────────────────────────────────
#' 使用 sf::st_distance 计算经纬度点之间的测地线距离（km）
#' 修复 v2 中 dist() 在经纬度上用欧氏距离的 bug (12:57-58)
#'
#' @param lon 经度向量
#' @param lat 纬度向量
#' @param unit 距离单位："km" 或 "m"
#' @return 距离矩阵 (n × n)，单位 km
geodesic_dist_matrix <- function(lon, lat, unit = "km") {
  stopifnot(length(lon) == length(lat))
  sf::sf_use_s2(FALSE)  # 使用 GEOS 而非 S2（更精确的测地线）

  pts <- sf::st_point_on_surface(
    sf::st_as_sf(
      data.frame(lon = lon, lat = lat),
      coords = c("lon", "lat"),
      crs = sf::st_crs(4326)
    )
  )

  d <- sf::st_distance(pts, pts)
  units(d) <- NULL  # 去除 units 属性，方便数值操作

  if (unit == "km") d <- d / 1000
  d
}

# ── 投影坐标到适合中国的 Albers 投影 ─────────────────────────────────
#' 中国 Albers 等面积投影（标准参数）
project_china_albers <- function(sf_obj) {
  sf::st_transform(sf_obj, crs = sf::st_crs(
    "+proj=aea +lat_1=25 +lat_2=47 +lat_0=30 +lon_0=105 +
     x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  ))
}

# ── 从经纬度生成 sf 点 ───────────────────────────────────────────────
make_sf_points <- function(lon, lat, crs = 4326) {
  sf::st_as_sf(
    data.frame(lon = lon, lat = lat),
    coords = c("lon", "lat"),
    crs = sf::st_crs(crs)
  )
}

# ── 计算网格质心间距离（快捷版）────────────────────────────────────────
#' 给定 grid_env data.frame（含 centroid_lon, centroid_lat），
#' 返回 km 距离矩阵
grid_distance_matrix <- function(grid_env) {
  geodesic_dist_matrix(grid_env$centroid_lon, grid_env$centroid_lat, unit = "km")
}

message("[utils_spatial] loaded")
