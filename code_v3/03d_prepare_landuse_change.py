#!/usr/bin/env python3
"""03d_prepare_landuse_change.py — v3 土地利用变化数据准备 (v2: 窗口裁剪法)

核心优化：对 135079×161378 的大型 CLCD GeoTIFF，不做全局 zonal_stats，
而是对每个 100km 网格：先用 rasterio.window 裁剪小窗口，再统计类别频率。
预计每个年份 1247 网格 × ~5s = ~2h。

输入:
  - CLCD GeoTIFF 文件 (CLCD_v01_YYYY_albert.tif, Albers投影, uint8)
  - 100km 网格 GeoJSON (china_grid_100km_active.geojson, WGS84)

输出:
  - data/derived_v3/landuse_change_v3.rds
  - data/derived_v3/landuse_by_year_v3.rds
  - results_v3/table_landuse_change_v3.csv
"""

import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import from_bounds
import geopandas as gpd

warnings.filterwarnings('ignore')

# ── 配置 ──────────────────────────────────────────────────────────
PROJECT_ROOT = os.environ.get(
    'V3_PROJECT_ROOT',
    os.path.expanduser('~/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3')
)

CLCD_DIR = os.environ.get(
    'V3_CLCD_DIR',
    os.path.join(PROJECT_ROOT, 'data', 'external', 'clcd')
)

CLCD_DIR_OLD = os.path.expanduser('~/bird_dynamic_occupancy_analysis/data/external/clcd')

DERIVED_DIR = os.path.join(PROJECT_ROOT, 'data', 'derived_v3')
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results_v3')
GRID_GEOJSON = os.path.join(DERIVED_DIR, 'china_grid_100km_active.geojson')

# 需要的关键年份 (P1-P5 对应 2000, 2005, 2010, 2015, 2020)
CLCD_YEARS = [2000, 2005, 2010, 2015, 2020]

# CLD 分类编码
CLCD_CLASSES = {
    'cropland': 1, 'forest': 2, 'shrub': 3, 'grassland': 4,
    'water': 5, 'snow_ice': 6, 'barren': 7, 'impervious': 8, 'wetland': 9
}

# 生态意义分组
NATURAL_CLASSES = ['forest', 'shrub', 'grassland', 'water', 'wetland']

print(f"[03d] Starting land use change data preparation (v2: windowed extraction)")
t0 = time.time()

# ── 1. 加载网格 ──────────────────────────────────────────────────
print("[03d] Loading grid...")
grid = gpd.read_file(GRID_GEOJSON)
print(f"[03d] Grid: {len(grid)} cells")

if 'centroid_lon' not in grid.columns:
    centroids = grid.geometry.centroid
    grid['centroid_lon'] = centroids.x
    grid['centroid_lat'] = centroids.y

# ── 2. 查找 CLCD 文件 ────────────────────────────────────────────
clcd_tif_paths = {}
for yr in CLCD_YEARS:
    candidates = [
        os.path.join(CLCD_DIR, f'CLCD_v01_{yr}_albert.tif'),
        os.path.join(CLCD_DIR_OLD, f'CLCD_v01_{yr}_albert.tif'),
        os.path.join(CLCD_DIR, f'clcd_{yr}.tif'),
    ]
    for c in candidates:
        real = os.path.realpath(c) if os.path.islink(c) else c
        if os.path.exists(real):
            clcd_tif_paths[yr] = c  # 保留原始路径(可能含软链接)
            break

available_years = sorted(clcd_tif_paths.keys())
print(f"[03d] Available CLCD years: {available_years}")

if len(available_years) < 2:
    print(f"[03d] ERROR: Need at least 2 CLCD years, found {len(available_years)}")
    sys.exit(1)

# ── 3. 窗口裁剪法提取各年份各网格的土地覆盖面积比例 ──────────────
def extract_landcover_windowed(tif_path, grid_gdf, class_codes):
    """
    对每个网格：裁剪对应窗口 → 读入小数组 → np.bincount 统计频率。
    比全局 zonal_stats 快 10-50 倍，因为只读必要的像素。
    """
    print(f"[03d] Reading {os.path.basename(tif_path)}")

    grid_proj = grid_gdf.to_crs('EPSG:4326')  # 先确保 WGS84
    # 读取栅格 CRS，将网格投影到栅格 CRS
    with rasterio.open(tif_path) as src:
        raster_crs = src.crs
        raster_shape = src.shape
        print(f"[03d]   CRS: {raster_crs}")
        print(f"[03d]   Shape: {raster_shape}")

    grid_in_raster_crs = grid_gdf.to_crs(raster_crs)

    n_cells = len(grid_in_raster_crs)
    # 预分配结果
    class_names = list(class_codes.keys())
    result = {cn: np.full(n_cells, np.nan) for cn in class_names}
    result['grid_cell'] = grid_gdf['grid_cell'].values

    t_start = time.time()
    n_done = 0

    with rasterio.open(tif_path) as src:
        for idx, row in grid_in_raster_crs.iterrows():
            geom = row.geometry
            minx, miny, maxx, maxy = geom.bounds

            # 加 1km 缓冲确保覆盖（100km 网格边缘像素）
            buf = 1000  # 1km in Albers meters
            minx -= buf; miny -= buf; maxx += buf; maxy += buf

            try:
                window = from_bounds(minx, miny, maxx, maxy,
                                     src.transform)
                # 只读窗口区域
                data = src.read(1, window=window)

                if data.size == 0:
                    continue

                # 忽略 NoData
                valid = data[data != src.nodata] if src.nodata is not None else data.ravel()
                if valid.size == 0:
                    continue

                total = valid.size
                counts = np.bincount(valid.ravel(), minlength=10)

                for cn, cv in class_codes.items():
                    if cv < len(counts):
                        result[cn][idx] = counts[cv] / total

            except Exception as e:
                # 某些边缘网格可能超出栅格范围，跳过
                pass

            n_done += 1
            if n_done % 100 == 0:
                elapsed = time.time() - t_start
                rate = n_done / elapsed
                eta = (n_cells - n_done) / rate if rate > 0 else 0
                print(f"[03d]   {n_done}/{n_cells} cells ({rate:.1f} cells/s, ETA {eta:.0f}s)")

    elapsed = time.time() - t_start
    print(f"[03d]   Extracted {n_cells} cells in {elapsed:.1f}s ({n_cells/elapsed:.1f} cells/s)")

    return pd.DataFrame(result)


# 提取各年份
landcover_by_year = {}
for yr in available_years:
    tif_path = clcd_tif_paths[yr]
    lc = extract_landcover_windowed(tif_path, grid, CLCD_CLASSES)

    # 重命名列加年份后缀
    rename_map = {c: f'{c}_{yr}' for c in CLCD_CLASSES.keys()}
    lc = lc.rename(columns=rename_map)

    landcover_by_year[yr] = lc
    print(f"[03d] Year {yr}: {len(lc)} grids extracted")

# 合并所有年份
lc_all = landcover_by_year[available_years[0]]
for yr in available_years[1:]:
    lc_all = lc_all.merge(landcover_by_year[yr], on='grid_cell', how='outer')

print(f"[03d] Combined landcover data: {lc_all.shape}")

# ── 4. 计算各时期变化量 (P1→P2, P2→P3, ..., 以及总变化 P1→P5) ──
early_year = min(available_years)
late_year = max(available_years)
print(f"[03d] Computing land use change: {early_year} -> {late_year}")

lc_types = ['forest', 'shrub', 'grassland', 'wetland', 'water', 'impervious', 'cropland']

landuse_change = pd.DataFrame({'grid_cell': lc_all['grid_cell']})

# 总变化 (early → late)
for lt in lc_types:
    early_col = f'{lt}_{early_year}'
    late_col = f'{lt}_{late_year}'
    if early_col in lc_all.columns and late_col in lc_all.columns:
        landuse_change[f'delta_{lt}'] = lc_all[late_col] - lc_all[early_col]
    else:
        landuse_change[f'delta_{lt}'] = np.nan

# 自然地类总比例变化
natural_early_cols = [f'{c}_{early_year}' for c in NATURAL_CLASSES if f'{c}_{early_year}' in lc_all.columns]
natural_late_cols = [f'{c}_{late_year}' for c in NATURAL_CLASSES if f'{c}_{late_year}' in lc_all.columns]
if natural_early_cols and natural_late_cols:
    natural_early = lc_all[natural_early_cols].sum(axis=1, skipna=True)
    natural_late = lc_all[natural_late_cols].sum(axis=1, skipna=True)
    landuse_change['delta_natural'] = natural_late - natural_early
else:
    landuse_change['delta_natural'] = np.nan

# 各时期比例 (P1-P5)
period_year_map = {
    'P1': 2000, 'P2': 2005, 'P3': 2010, 'P4': 2015, 'P5': 2020
}
for period, yr in period_year_map.items():
    if yr in available_years:
        for lt in lc_types:
            col = f'{lt}_{yr}'
            if col in lc_all.columns:
                landuse_change[f'{lt}_{period}'] = lc_all[col]

# 添加坐标和元数据
grid_idx = grid.set_index('grid_cell')
landuse_change['centroid_lon'] = grid_idx.loc[landuse_change['grid_cell'], 'centroid_lon'].values
landuse_change['centroid_lat'] = grid_idx.loc[landuse_change['grid_cell'], 'centroid_lat'].values
landuse_change['early_year'] = early_year
landuse_change['late_year'] = late_year
landuse_change['data_source'] = 'CLCD'

# ── 5. 准备 landuse_full (含各年份面积比例) ─────────────────────
landuse_full = lc_all.copy()
landuse_full['centroid_lon'] = grid_idx.loc[landuse_full['grid_cell'], 'centroid_lon'].values
landuse_full['centroid_lat'] = grid_idx.loc[landuse_full['grid_cell'], 'centroid_lat'].values
cols_full = ['grid_cell', 'centroid_lon', 'centroid_lat'] + \
    [c for c in landuse_full.columns if c not in ['grid_cell', 'centroid_lon', 'centroid_lat']]
landuse_full = landuse_full[cols_full]

# ── 6. 保存 ──────────────────────────────────────────────────────
os.makedirs(RESULTS_DIR, exist_ok=True)

csv_path = os.path.join(DERIVED_DIR, 'landuse_change_v3.csv')
landuse_change.to_csv(csv_path, index=False)
print(f"[03d] Saved CSV: {csv_path}")

csv_results = os.path.join(RESULTS_DIR, 'table_landuse_change_v3.csv')
landuse_change.to_csv(csv_results, index=False)
print(f"[03d] Saved results CSV: {csv_results}")

# 保存 RDS
try:
    import pyreadr
    pyreadr.write_rds(os.path.join(DERIVED_DIR, 'landuse_change_v3.rds'), landuse_change)
    pyreadr.write_rds(os.path.join(DERIVED_DIR, 'landuse_by_year_v3.rds'), landuse_full)
    print(f"[03d] Saved RDS files")
except Exception as e:
    print(f"[03d] pyreadr RDS save failed: {e}")

# ── 汇总 ─────────────────────────────────────────────────────────
print(f"[03d] {len(landuse_change)} grid cells with landuse change data")
print(f"[03d] Period: {early_year} -> {late_year}")
for lt in lc_types:
    delta_col = f'delta_{lt}'
    if delta_col in landuse_change.columns:
        valid = landuse_change[delta_col].dropna()
        if len(valid) > 0:
            print(f"[03d]   {lt}: [{valid.min():.4f}, {valid.max():.4f}] (mean {valid.mean():.4f})")

elapsed = time.time() - t0
print(f"[03d] Completed in {elapsed:.1f}s")
