#!/usr/bin/env python3
"""03c_prepare_climate_change.py — v3 气候变化数据准备

使用 Python netCDF4/xarray 读取 CRU TS NetCDF（terra GDAL 无法读取），
按 100km 网格聚合，计算 period 均值和气候变化量。

输入:
  - CRU TS v4.09 月值 NetCDF (tmp, tmx)
  - 100km 网格 GeoJSON (china_grid_100km_active.geojson)

输出:
  - data/derived_v3/climate_change_v3.rds  (via R conversion)
  - results_v3/table_climate_change_v3.csv
  - data/derived_v3/climate_change_v3.csv  (intermediate for R)
"""

import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
from rasterstats import zonal_stats

warnings.filterwarnings('ignore')

# ── 配置 ──────────────────────────────────────────────────────────
PROJECT_ROOT = os.environ.get(
    'V3_PROJECT_ROOT',
    os.path.expanduser('~/projects/bird-new-distribution-records/tasks/bird_dynamic_occupancy_analysis_v3')
)

CRU_DIR = os.environ.get(
    'V3_CRU_DIR',
    os.path.join(PROJECT_ROOT, 'data', 'external', 'cru_ts')
)

DERIVED_DIR = os.path.join(PROJECT_ROOT, 'data', 'derived_v3')
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results_v3')
GRID_GEOJSON = os.path.join(DERIVED_DIR, 'china_grid_100km_active.geojson')

# 研究期
ANALYSIS_YR_LO = 2000
ANALYSIS_YR_HI = 2024

# Period 定义
PERIODS = {
    'P1': (2000, 2004),
    'P2': (2005, 2009),
    'P3': (2010, 2014),
    'P4': (2015, 2019),
    'P5': (2020, 2024),
}

print(f"[03c] Starting climate change data preparation")
print(f"[03c] Project root: {PROJECT_ROOT}")
print(f"[03c] CRU dir: {CRU_DIR}")
t0 = time.time()

# ── 1. 加载网格 ──────────────────────────────────────────────────
print("[03c] Loading grid...")
grid = gpd.read_file(GRID_GEOJSON)
print(f"[03c] Grid: {len(grid)} cells")
# 确保经纬度列存在
if 'centroid_lon' not in grid.columns:
    centroids = grid.geometry.centroid
    grid['centroid_lon'] = centroids.x
    grid['centroid_lat'] = centroids.y

# ── 2. 查找 CRU TS 文件 ──────────────────────────────────────────
cru_files = {}
for varname in ['tmp', 'tmx']:
    # 尝试多种文件名模式
    candidates = [
        os.path.join(CRU_DIR, f'cru_ts4.09.1901.2024.{varname}.dat.nc'),
        os.path.join(CRU_DIR, f'cru_ts4.09.1901.2023.{varname}.dat.nc'),
    ]
    found = None
    for c in candidates:
        if os.path.exists(c):
            found = c
            break
        # 检查是否是软链接指向的文件存在
        if os.path.islink(c):
            real = os.path.realpath(c)
            if os.path.exists(real):
                found = real
                break

    if found:
        cru_files[varname] = found
        print(f"[03c] {varname}: {found}")
    else:
        print(f"[03c] WARNING: CRU TS {varname} not found in {CRU_DIR}")

if len(cru_files) == 0:
    print("[03c] ERROR: No CRU TS NetCDF files found!")
    sys.exit(1)

# ── 3. 提取各 period 均值 ────────────────────────────────────────
def extract_cru_period(nc_path, varname, grid_gdf, periods):
    """从 CRU TS NetCDF 提取各 period 的网格均值"""
    print(f"[03c] Reading CRU TS {varname} from {nc_path}")

    # 使用 xarray 打开 NetCDF（比 netCDF4 更方便处理时间维度）
    # CRU TS 使用 CF conventions，xarray 可以自动解析时间
    try:
        ds = xr.open_dataset(nc_path, decode_times=True)
    except Exception as e:
        print(f"[03c] xarray decode_times failed: {e}, trying without...")
        ds = xr.open_dataset(nc_path, decode_times=False)
        # 手动解析时间
        if 'time' in ds.dims:
            time_var = ds['time']
            if hasattr(time_var, 'units') and 'since' in str(time_var.units):
                import cftime
                times = cftime.num2date(time_var.values, time_var.units, calendar=time_var.attrs.get('calendar', 'standard'))
                ds['time'] = pd.DatetimeIndex([pd.Timestamp(t.year, t.month, t.day) for t in times])

    print(f"[03c] Variables: {list(ds.data_vars)}")
    print(f"[03c] Dims: {dict(ds.dims)}")

    # 获取数据变量名（通常与 varname 相同，但有时不同）
    data_var = varname
    if data_var not in ds.data_vars:
        # 尝试查找匹配的变量
        for v in ds.data_vars:
            if varname in v.lower():
                data_var = v
                break
    print(f"[03c] Using variable: {data_var}")

    # 时间范围
    if 'time' in ds.dims:
        time_vals = ds['time'].values
        year_min = pd.Timestamp(time_vals[0]).year
        year_max = pd.Timestamp(time_vals[-1]).year
        print(f"[03c] Time range: {year_min}-{year_max}")
    else:
        print(f"[03c] No time dimension found! Available dims: {list(ds.dims)}")
        ds.close()
        return None

    period_means = {}

    for pn, (yr_lo, yr_hi) in periods.items():
        yr_hi_actual = min(yr_hi, year_max, ANALYSIS_YR_HI)
        yr_lo_actual = max(yr_lo, year_min, ANALYSIS_YR_LO)

        if yr_lo_actual > yr_hi_actual:
            print(f"[03c] No data for {pn} ({yr_lo}-{yr_hi})")
            continue

        # 选择时间范围
        da = ds[data_var].sel(time=slice(f'{yr_lo_actual}-01-01', f'{yr_hi_actual}-12-31'))

        if len(da.time) == 0:
            print(f"[03c] No data for {pn}")
            continue

        # 计算年均值，然后 period 均值
        annual_mean = da.groupby('time.year').mean(dim='time')
        period_mean = annual_mean.mean(dim='year')

        print(f"[03c] {pn} ({yr_lo_actual}-{yr_hi_actual}): {len(da.time)} months, extracting values...")

        # CRU TS 是 0.5° 分辨率，100km 网格约覆盖 1° → 4 个 CRU 像素
        # 使用网格质心进行最近邻提取即可，简单高效
        lats = grid_gdf['centroid_lat'].values
        lons = grid_gdf['centroid_lon'].values

        # 批量提取：用 xarray 的 sel + 向量化
        try:
            # 方法1：用 xarray 的 interpolation（向量化，更快）
            period_mean_interp = period_mean.interp(
                lat=xr.DataArray(lats, dims='z'),
                lon=xr.DataArray(lons, dims='z'),
                method='nearest'
            )
            vals = period_mean_interp.values.astype(float)
            vals[np.isnan(vals)] = np.nan
        except Exception as e:
            # 方法2：逐点提取（更可靠但较慢）
            print(f"[03c] Vectorized extraction failed ({e}), using point-by-point...")
            vals = []
            for lon, lat in zip(lons, lats):
                try:
                    v = float(period_mean.sel(lat=lat, lon=lon, method='nearest').values)
                    vals.append(v if not np.isnan(v) else np.nan)
                except Exception:
                    vals.append(np.nan)

        period_means[pn] = list(vals)
        mean_val = np.nanmean(vals)
        n_valid = np.sum(~np.isnan(np.array(vals, dtype=float)))
        print(f"[03c] {pn}: mean = {mean_val:.2f} ({n_valid}/{len(vals)} valid)")

    ds.close()
    return period_means


# 提取 mean temperature (tmp)
tmp_periods = None
if 'tmp' in cru_files:
    tmp_periods = extract_cru_period(cru_files['tmp'], 'tmp', grid, PERIODS)

# 提取 max temperature (tmx)
tmx_periods = None
if 'tmx' in cru_files:
    tmx_periods = extract_cru_period(cru_files['tmx'], 'tmx', grid, PERIODS)

# ── 4. 计算气候变化量 ────────────────────────────────────────────
if tmp_periods is None and tmx_periods is None:
    print("[03c] ERROR: No climate data extracted!")
    sys.exit(1)

climate_data = {'grid_cell': grid['grid_cell'].values}

# 处理 tmp 数据
if tmp_periods:
    # 历史基线（P1+P2 均值）
    baseline_periods = [p for p in ['P1', 'P2'] if p in tmp_periods]
    if baseline_periods:
        tmp_baseline = np.nanmean([tmp_periods[p] for p in baseline_periods], axis=0)
        tmp_baseline_sd = np.nanstd([tmp_periods[p] for p in baseline_periods], axis=0)
    else:
        # 使用最早可用 period
        first_p = list(tmp_periods.keys())[0]
        tmp_baseline = np.array(tmp_periods[first_p])
        tmp_baseline_sd = np.ones(len(tmp_baseline)) * 0.5  # 默认 SD

    # 当前值（最晚可用 period）
    current_periods = [p for p in ['P5', 'P4', 'P3'] if p in tmp_periods]
    tmp_current = np.array(tmp_periods[current_periods[0]])

    climate_data['delta_t_mean'] = tmp_current - tmp_baseline
    climate_data['delta_t_std'] = (tmp_current - tmp_baseline) / np.maximum(tmp_baseline_sd, 0.01)

    # 保存各 period 均值
    for pn in ['P1', 'P2', 'P3', 'P4', 'P5']:
        key = f'tmp_{pn}'
        if pn in tmp_periods:
            climate_data[key] = tmp_periods[pn]
        else:
            climate_data[key] = np.nan
else:
    climate_data['delta_t_mean'] = np.nan
    climate_data['delta_t_std'] = np.nan
    for pn in ['P1', 'P2', 'P3', 'P4', 'P5']:
        climate_data[f'tmp_{pn}'] = np.nan

# 处理 tmx 数据
if tmx_periods:
    baseline_periods = [p for p in ['P1', 'P2'] if p in tmx_periods]
    if baseline_periods:
        tmx_baseline = np.nanmean([tmx_periods[p] for p in baseline_periods], axis=0)
    else:
        first_p = list(tmx_periods.keys())[0]
        tmx_baseline = np.array(tmx_periods[first_p])

    current_periods = [p for p in ['P5', 'P4', 'P3'] if p in tmx_periods]
    tmx_current = np.array(tmx_periods[current_periods[0]])

    climate_data['delta_t_extreme'] = tmx_current - tmx_baseline

    for pn in ['P1', 'P2', 'P3', 'P4', 'P5']:
        key = f'tmx_{pn}'
        if pn in tmx_periods:
            climate_data[key] = tmx_periods[pn]
        else:
            climate_data[key] = np.nan
else:
    climate_data['delta_t_extreme'] = np.nan
    for pn in ['P1', 'P2', 'P3', 'P4', 'P5']:
        climate_data[f'tmx_{pn}'] = np.nan

# 添加坐标
climate_data['centroid_lon'] = grid['centroid_lon'].values
climate_data['centroid_lat'] = grid['centroid_lat'].values

# 创建 DataFrame
climate_df = pd.DataFrame(climate_data)
# 重新排列列
cols = ['grid_cell', 'centroid_lon', 'centroid_lat',
        'delta_t_mean', 'delta_t_extreme', 'delta_t_std',
        'tmp_P1', 'tmp_P2', 'tmp_P3', 'tmp_P4', 'tmp_P5',
        'tmx_P1', 'tmx_P2', 'tmx_P3', 'tmx_P4', 'tmx_P5']
climate_df = climate_df[[c for c in cols if c in climate_df.columns]]

# ── 5. 保存 ──────────────────────────────────────────────────────
os.makedirs(RESULTS_DIR, exist_ok=True)

csv_path = os.path.join(DERIVED_DIR, 'climate_change_v3.csv')
climate_df.to_csv(csv_path, index=False)
print(f"[03c] Saved CSV: {csv_path}")

csv_results = os.path.join(RESULTS_DIR, 'table_climate_change_v3.csv')
climate_df.to_csv(csv_results, index=False)
print(f"[03c] Saved results CSV: {csv_results}")

# 保存为 RDS（通过 pyreadr）
try:
    import pyreadr
    rds_path = os.path.join(DERIVED_DIR, 'climate_change_v3.rds')
    pyreadr.write_rds(rds_path, climate_df)
    print(f"[03c] Saved RDS: {rds_path}")
except Exception as e:
    print(f"[03c] pyreadr RDS save failed: {e}")
    print(f"[03c] CSV saved, use R to convert to RDS")

# ── 汇总 ─────────────────────────────────────────────────────────
print(f"[03c] {len(climate_df)} grid cells with climate change data")
if 'delta_t_mean' in climate_df.columns:
    valid = climate_df['delta_t_mean'].dropna()
    if len(valid) > 0:
        print(f"[03c] delta_t_mean range: [{valid.min():.2f}, {valid.max():.2f}]")
if 'delta_t_extreme' in climate_df.columns:
    valid = climate_df['delta_t_extreme'].dropna()
    if len(valid) > 0:
        print(f"[03c] delta_t_extreme range: [{valid.min():.2f}, {valid.max():.2f}]")

elapsed = time.time() - t0
print(f"[03c] Completed in {elapsed:.1f}s")
