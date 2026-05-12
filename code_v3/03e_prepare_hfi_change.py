#!/usr/bin/env python3
"""03e_prepare_hfi_change.py — v3 人类足迹指数变化数据准备

从本地 HFI 栅格 (data/external/hfi/) 提取各 100km 网格的 HFI 值，
计算 delta_hfi (2020 - 2000) 及各 period 的 HFI 值。

数据源: Mu et al. 2022, Scientific Data
  - 年度全球 HFI, ~1km 分辨率, Mollweide 投影
  - 本地已有: hfp2000.zip, hfp2005.zip, ..., hfp2024.zip

输出:
  - data/derived_v3/hfi_change_v3.rds
  - data/derived_v3/hfi_by_year_v3.rds
  - results_v3/table_hfi_change_v3.csv
"""

import os
import sys
import time
import warnings
import zipfile
import tempfile
import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask as rio_mask
from rasterstats import zonal_stats
import geopandas as gpd
from shapely.geometry import box

warnings.filterwarnings('ignore')

# ── 配置 ──────────────────────────────────────────────────────────
PROJECT_ROOT = os.environ.get(
    'V3_PROJECT_ROOT',
    os.path.expanduser('~/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3')
)

HFI_DIR = os.path.join(PROJECT_ROOT, 'data', 'external', 'hfi')
DERIVED_DIR = os.path.join(PROJECT_ROOT, 'data', 'derived_v3')
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results_v3')
GRID_GEOJSON = os.path.join(DERIVED_DIR, 'china_grid_100km_active.geojson')

# 需要 HFI 数据的关键年份
HFI_YEARS = [2000, 2005, 2010, 2015, 2020, 2024]

# Period 映射
PERIOD_MAP = {
    2000: 'P1', 2005: 'P2', 2010: 'P3', 2015: 'P4', 2020: 'P5', 2024: 'P5'
}

print(f"[03e] Starting HFI change data preparation")
t0 = time.time()

# ── 1. 加载网格 ──────────────────────────────────────────────────
print("[03e] Loading grid...")
grid = gpd.read_file(GRID_GEOJSON)
print(f"[03e] Grid: {len(grid)} cells")

if 'centroid_lon' not in grid.columns:
    centroids = grid.geometry.centroid
    grid['centroid_lon'] = centroids.x
    grid['centroid_lat'] = centroids.y

# ── 2. 查找/解压 HFI 文件 ───────────────────────────────────────
hfi_tif_paths = {}

for yr in HFI_YEARS:
    # 先检查是否已有解压的 tif
    tif_path = os.path.join(HFI_DIR, f'hfp{yr}.tif')
    if os.path.exists(tif_path):
        hfi_tif_paths[yr] = tif_path
        print(f"[03e] {yr}: already extracted {tif_path}")
        continue

    # 检查 zip 文件
    zip_path = os.path.join(HFI_DIR, f'hfp{yr}.zip')
    if os.path.exists(zip_path):
        print(f"[03e] {yr}: extracting from {zip_path}...")
        try:
            with zipfile.ZipFile(zip_path, 'r') as zf:
                # 获取 zip 内的 tif 文件名
                tif_names = [n for n in zf.namelist() if n.endswith('.tif')]
                if tif_names:
                    zf.extract(tif_names[0], HFI_DIR)
                    extracted_path = os.path.join(HFI_DIR, tif_names[0])
                    # 如果解压路径与预期不同，重命名
                    if extracted_path != tif_path:
                        os.rename(extracted_path, tif_path)
                    hfi_tif_paths[yr] = tif_path
                    print(f"[03e] {yr}: extracted to {tif_path}")
                else:
                    print(f"[03e] {yr}: no .tif found in zip!")
        except Exception as e:
            print(f"[03e] {yr}: extraction failed: {e}")
    else:
        print(f"[03e] {yr}: zip file not found at {zip_path}")

available_years = sorted(hfi_tif_paths.keys())
print(f"[03e] Available HFI years: {available_years}")

if len(available_years) < 2:
    print(f"[03e] ERROR: Need at least 2 HFI years, found {len(available_years)}")
    sys.exit(1)

# ── 3. 从栅格提取各年份 HFI ──────────────────────────────────────
def extract_hfi_from_raster(tif_path, grid_gdf):
    """从 HFI 栅格提取各网格的 HFI 均值"""
    file_size_mb = os.path.getsize(tif_path) / 1e6
    print(f"[03e]     Reading {os.path.basename(tif_path)} ({file_size_mb:.0f} MB)...")

    with rasterio.open(tif_path) as src:
        raster_crs = src.crs
        print(f"[03e]     CRS: {raster_crs}")

    # 将网格投影到栅格 CRS (Mollweide)
    grid_proj = grid_gdf.to_crs(raster_crs)

    # 裁剪到网格范围以减少内存占用
    # 获取网格总范围的 bbox
    total_bounds = grid_proj.total_bounds  # [minx, miny, maxx, maxy]

    # 使用 rasterstats 进行 zonal 统计
    print(f"[03e]     Computing zonal stats for {len(grid_proj)} cells...")
    stats = zonal_stats(
        grid_proj, tif_path,
        stats=['mean'],
        all_touched=True,
        nodata=0,  # HFI 中 0 表示无数据
    )

    vals = [s['mean'] if s['mean'] is not None else np.nan for s in stats]
    return vals


# 按文件大小排序（小文件先处理，快速验证流程）
year_sizes = {}
for yr in available_years:
    if yr in hfi_tif_paths:
        year_sizes[yr] = os.path.getsize(hfi_tif_paths[yr])
yr_order = sorted(available_years, key=lambda y: year_sizes.get(y, float('inf')))

hfi_data = {'grid_cell': grid['grid_cell'].values}

for yr in yr_order:
    tif_path = hfi_tif_paths[yr]
    if not os.path.exists(tif_path):
        continue

    vals = extract_hfi_from_raster(tif_path, grid)
    hfi_data[f'hfi_{yr}'] = vals

    n_valid = sum(1 for v in vals if not np.isnan(v))
    print(f"[03e]   Year {yr}: {n_valid}/{len(vals)} valid ({100*n_valid/len(vals):.1f}%)")

    # 释放内存
    import gc
    gc.collect()

# 创建 HFI DataFrame
hfi_all = pd.DataFrame(hfi_data)

# ── 4. 更新 grid_environment_v3.rds 中的 HFI 数据 ────────────────
grid_env_path = os.path.join(DERIVED_DIR, 'grid_environment_v3.rds')
if os.path.exists(grid_env_path):
    print(f"[03e] Updating grid_environment_v3.rds with new HFI values...")
    try:
        import pyreadr
        grid_env = pyreadr.read_rds(grid_env_path)
        if isinstance(grid_env, dict):
            grid_env = grid_env[list(grid_env.keys())[0]]

        # 移除旧 HFI 列
        old_hfi_cols = [c for c in grid_env.columns if c.startswith('hfi_')]
        grid_env = grid_env.drop(columns=old_hfi_cols, errors='ignore')

        # 合并新 HFI 数据
        grid_env = grid_env.merge(hfi_all, on='grid_cell', how='left')

        # 重新计算 hfi_mean 和 hfi_sd
        hfi_value_cols = [c for c in hfi_all.columns if c.startswith('hfi_') and c != 'grid_cell']
        if hfi_value_cols:
            grid_env['hfi_mean'] = grid_env[hfi_value_cols].mean(axis=1, skipna=True)
            grid_env['hfi_sd'] = grid_env[hfi_value_cols].std(axis=1, skipna=True)

        pyreadr.write_rds(grid_env_path, grid_env)
        n_valid = grid_env['hfi_2000'].notna().sum() if 'hfi_2000' in grid_env.columns else 0
        print(f"[03e] grid_environment_v3 updated: hfi_2000 valid = {n_valid}")
    except Exception as e:
        print(f"[03e] grid_environment update failed: {e}")

# ── 5. 计算 HFI 变化量 ──────────────────────────────────────────
available_hfi_cols = [c for c in hfi_all.columns if c.startswith('hfi_') and c != 'grid_cell']
available_hfi_years = sorted([int(c.replace('hfi_', '')) for c in available_hfi_cols])

early_yr = min(available_hfi_years)
late_yr = max(available_hfi_years)

early_col = f'hfi_{early_yr}'
late_col = f'hfi_{late_yr}'

print(f"[03e] Computing HFI change: {early_yr} -> {late_yr}")

hfi_change = hfi_all.copy()
hfi_change['delta_hfi'] = hfi_change[late_col] - hfi_change[early_col]
hfi_change['early_year'] = early_yr
hfi_change['late_year'] = late_yr
hfi_change['data_source'] = 'HFI_Mu_2022'

# 添加坐标
hfi_change = hfi_change.merge(
    grid[['grid_cell', 'centroid_lon', 'centroid_lat']],
    on='grid_cell', how='left'
)

# ── 6. 计算各 period 的 HFI 均值 ────────────────────────────────
period_hfi_long = hfi_all.melt(
    id_vars='grid_cell',
    value_vars=available_hfi_cols,
    var_name='year_str',
    value_name='hfi'
)
period_hfi_long['year'] = period_hfi_long['year_str'].str.replace('hfi_', '').astype(int)
period_hfi_long['period'] = period_hfi_long['year'].map(PERIOD_MAP)
period_hfi_long = period_hfi_long.dropna(subset=['period'])

period_hfi = period_hfi_long.groupby(['grid_cell', 'period'])['hfi'].mean().reset_index()
period_hfi = period_hfi.pivot(index='grid_cell', columns='period', values='hfi')
period_hfi.columns = [f'hfi_{c}' for c in period_hfi.columns]
period_hfi = period_hfi.reset_index()

hfi_change = hfi_change.merge(period_hfi, on='grid_cell', how='left')

# 重新排列列
cols = ['grid_cell', 'centroid_lon', 'centroid_lat'] + \
       [c for c in hfi_change.columns if c.startswith('hfi_') and c != 'grid_cell'] + \
       ['delta_hfi', 'early_year', 'late_year', 'data_source']
cols = [c for c in cols if c in hfi_change.columns]
hfi_change = hfi_change[cols]

# ── 7. 保存 ──────────────────────────────────────────────────────
os.makedirs(RESULTS_DIR, exist_ok=True)

csv_path = os.path.join(DERIVED_DIR, 'hfi_change_v3.csv')
hfi_change.to_csv(csv_path, index=False)
print(f"[03e] Saved CSV: {csv_path}")

csv_results = os.path.join(RESULTS_DIR, 'table_hfi_change_v3.csv')
hfi_change.to_csv(csv_results, index=False)
print(f"[03e] Saved results CSV: {csv_results}")

# 保存 RDS
try:
    import pyreadr
    pyreadr.write_rds(os.path.join(DERIVED_DIR, 'hfi_change_v3.rds'), hfi_change)
    pyreadr.write_rds(os.path.join(DERIVED_DIR, 'hfi_by_year_v3.rds'), hfi_all)
    print(f"[03e] Saved RDS files")
except Exception as e:
    print(f"[03e] pyreadr RDS save failed: {e}")

# ── 汇总 ─────────────────────────────────────────────────────────
print(f"[03e] {len(hfi_change)} grid cells with HFI change data")
print(f"[03e] Period: {early_yr} -> {late_yr}")
if 'delta_hfi' in hfi_change.columns:
    valid = hfi_change['delta_hfi'].dropna()
    if len(valid) > 0:
        print(f"[03e] delta_hfi range: [{valid.min():.2f}, {valid.max():.2f}] (mean {valid.mean():.2f})")
        print(f"[03e] N valid delta_hfi: {len(valid)} / {len(hfi_change)}")

elapsed = time.time() - t0
print(f"[03e] Completed in {elapsed:.1f}s")
