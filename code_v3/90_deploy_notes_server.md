# v3 服务器部署说明 / Server deployment notes

## 目标 / Goal
将 `bird_dynamic_occupancy_analysis/code_v3` 与其最小必要依赖同步到服务器，并在服务器端运行完整的 `v3` 动态占域分析工作流。

## 最小同步集合 / Minimal sync set
以下文件或目录已经整理在 [90_sync_manifest_v3.txt](/Users/dingchenchen/Documents/New%20project/bird_dynamic_occupancy_analysis/code_v3/90_sync_manifest_v3.txt)：

- `code_v3`
- `results_v2`
- `data/china_boundary`
- `data/中国shp`
- `data/derived_v2/china_grid_100km_v2.rds`
- `data/derived_v2/combined_events_merged_dedup_2000_2025.rds`
- `data/derived_v2/grid_environment_dynamic_occupancy.rds`
- `data/derived_v2/species_visit_2000_2024.rds`
- `data/derived_v2/visit_effort_2000_2024.rds`
- `data/external/cru_ts`
- `data/external/clcd`
- `data/external/hfi`

## 说明 / Notes
- 这份清单避免了把 `data/derived_v2` 里不再需要的超大旧模型对象一并上传。
- `AVONET` 与 `EltonTraits` 默认仍按 `00_config.R` 的服务器环境变量/默认路径读取；若服务器缺失，需要在服务器上补充对应路径或设置 `V3_AVONET_PATH`、`V3_ELTONTRAITS_PATH`。
- 服务器运行入口脚本为 [90_server_run_v3_pipeline.sh](/Users/dingchenchen/Documents/New%20project/bird_dynamic_occupancy_analysis/code_v3/90_server_run_v3_pipeline.sh)。
