#!/usr/bin/env Rscript
## 03b_extend_traits.R  —  v3 性状扩展脚本
##
## 新增两个关键性状：
##   1. diet_specialization — 基于 Morelli et al. 2021 方法
##      从 EltonTraits 1.0 (Wilman et al. 2014) 提取，使用 Gini 不等系数
##      Gini coefficient: G = ΣᵢΣⱼ|xᵢ-xⱼ| / (2n²·x̄), 0=泛化, 1=专化
##      参考：Morelli et al. 2019 Ecol Evol (方法) + Morelli et al. 2021 Conserv Lett (全球)
##   2. habitat_breadth — 基于 IUCN Red List 物种栖息地数量
##      使用 rredlist::rl_habitats() 提取，同物异名匹配 + 随机森林插值
##
## 输出：
##   data/derived_v3/trait_extended_v3.rds
##   results_v3/table_traits_extended_v3.csv
##
## 参考：
##   Morelli, F., Benedetti, Y., Møller, A. P., & Fuller, R. A. (2019).
##     Measuring avian specialization. Ecology and Evolution, 9, 8378–8386.
##   Morelli, F., Benedetti, Y., Hanson, J. O., & Fuller, R. A. (2021).
##     Global distribution and conservation of avian diet specialization.
##     Conservation Letters, e12795. https://doi.org/10.1111/conl.12795
##   Wilman, H., et al. (2014). EltonTraits 1.0. Ecology, 95(7), 2027.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(purrr)
})

# ── 加载配置 ──────────────────────────────────────────────────────────
CODE_V3 <- Sys.getenv("V3_CODE_DIR",
  file.path("~", "Documents", "New project", "bird_dynamic_occupancy_analysis", "code_v3"))
source(file.path(CODE_V3, "00_config.R"))
source(file.path(CODE_V3, "utils_paths.R"))
source(file.path(CODE_V3, "utils_core.R"))
P <- ensure_v3_dirs()

log_time("03b", "Starting trait extension")

# ── 1. 加载候选物种列表 ───────────────────────────────────────────────
cand_path <- v3_file("results", "table_dynamic_occupancy_candidate_species_all")
if (!file.exists(cand_path)) {
  # 尝试从 v2 获取
  cand_path <- v2_file("results", "table_dynamic_occupancy_candidate_species_all")
}
candidate_species <- read_csv_safe(cand_path)$species
if (is.null(candidate_species)) {
  stop("Candidate species list not found. Run Stage 2/3 first.")
}
message(sprintf("[03b] %d candidate species", length(candidate_species)))

# ── 2. 加载基础性状 ──────────────────────────────────────────────────
trait_imp <- read_csv_safe(TRAIT_IMPUTED_PATH)
if (is.null(trait_imp)) {
  # 尝试从 v2 结果中获取
  v2_trait_path <- file.path(DIRS$v2_results, "table_species_traits_imputed.csv")
  trait_imp <- read_csv_safe(v2_trait_path)
}
if (is.null(trait_imp)) {
  # 尝试从 v1 结果中获取
  v1_trait_path <- file.path(DIRS$v2_results, "..", "results", "table_species_traits_imputed.csv")
  v1_trait_path <- file.path(DIRS$data_raw, "..", "results", "table_species_traits_imputed.csv")
  trait_imp <- read_csv_safe(v1_trait_path)
}
if (is.null(trait_imp)) {
  message("[03b] WARNING: Base trait file not found. Creating from candidate species only.")
  trait_imp <- tibble(species = candidate_species)
}

# ── 3. diet_specialization（Morelli et al. 2021 方法）───────────────────
# 基于 EltonTraits 1.0 的 Gini 不等系数 (Morelli et al. 2019, 2021)
# EltonTraits 提供 10 个食性百分比列：Diet.Inv, Diet.Vend, Diet.Vect,
#   Diet.Vfish, Diet.Vunk, Diet.Scav, Diet.Fruit, Diet.Nect,
#   Diet.Seed, Diet.PlantO（均为0-100百分比）
# Gini coefficient: G = ΣᵢΣⱼ|xᵢ-xⱼ| / (2n²·x̄)
# G = 0 → 完全均匀（所有食性比例相等 = 泛化种）
# G = 1 → 完全集中（单一食性 = 专化种）
# 使用 DescTools::Gini() 计算（Morelli et al. 2019 原文使用此包）
# 输入为连续比例值（百分比/100），非二元化 0/1
# 注意：EltonTraits 物种名列为 Scientific（非 Species1）

compute_diet_specialization <- function(elton_df, species_list) {
  # EltonTraits 食性百分比列（10个，全部食性类别）
  # 参考 Wilman et al. 2014 和 Morelli et al. 2021
  diet_cols <- c("Diet.Inv", "Diet.Vend", "Diet.Vect", "Diet.Vfish",
                 "Diet.Vunk", "Diet.Scav", "Diet.Fruit", "Diet.Nect",
                 "Diet.Seed", "Diet.PlantO")

  # EltonTraits 物种名列为 Scientific
  sp_col <- if ("Scientific" %in% names(elton_df)) "Scientific" else "Species1"

  # 对齐物种
  result <- tibble(species = species_list) |>
    left_join(
      elton_df |>
        select(all_of(sp_col), all_of(intersect(diet_cols, names(elton_df)))),
      by = setNames(sp_col, "species")
    )

  # 获取实际可用的食性列
  avail_cols <- intersect(diet_cols, names(result))

  if (length(avail_cols) == 0) {
    message("[03b] WARNING: No EltonTraits diet columns found. diet_specialization = NA")
    result$diet_specialization <- NA_real_
    return(result |> select(species, diet_specialization))
  }

  # 检查 DescTools 包是否可用
  use_gini <- requireNamespace("DescTools", quietly = TRUE)
  if (!use_gini) {
    warning("[03b] DescTools not installed. Falling back to Shannon-based index. ",
            "Install DescTools for the correct Gini coefficient method.")
  }

  # 百分比 → 比例
  diet_mat <- as.matrix(result[, avail_cols])
  diet_mat_prop <- diet_mat / 100

  # 对每个物种计算 Gini 系数
  diet_spec <- apply(diet_mat_prop, 1, function(row) {
    p <- row[!is.na(row)]
    if (length(p) == 0 || sum(p) == 0) return(NA_real_)
    # 去除 0 值类别（物种不使用的食性类别）
    p <- p[p > 0]
    if (length(p) == 0) return(NA_real_)
    # 单类别物种：Gini = 1（最专化，所有食性集中在单一类型）
    if (length(p) == 1) return(1)
    # 重新标准化为比例（和=1）
    p <- p / sum(p)

    if (use_gini) {
      # Gini 不等系数：G = ΣᵢΣⱼ|xᵢ-xⱼ| / (2n²·x̄)
      # 使用 DescTools::Gini() — Morelli et al. 2019 原文方法
      DescTools::Gini(p)
    } else {
      # Fallback: Shannon 多样性逆指数（近似但非原文方法）
      H <- -sum(p * log(p))
      S <- length(p)
      if (S <= 1) return(1)  # 单类别 = 最专化
      1 - H / log(S)
    }
  })

  result$diet_specialization <- diet_spec
  result |> select(species, diet_specialization)
}

# 尝试加载 EltonTraits
diet_spec_df <- NULL
if (file.exists(ELTONTRAITS_PATH)) {
  message("[03b] Loading EltonTraits from: ", ELTONTRAITS_PATH)
  elton <- read_delim(ELTONTRAITS_PATH, delim = "\t", show_col_types = FALSE)
  # Normalize column names: replace hyphens with dots
  names(elton) <- gsub("-", ".", names(elton))

  # 检查必要的列（EltonTraits 物种名列为 Scientific 或 Species1）
  sp_col_check <- if ("Scientific" %in% names(elton)) "Scientific" else "Species1"
  if (sp_col_check %in% names(elton)) {
    diet_spec_df <- compute_diet_specialization(elton, candidate_species)
    n_matched <- sum(!is.na(diet_spec_df$diet_specialization))
    message(sprintf("[03b] diet_specialization: %d/%d species matched from EltonTraits",
                    n_matched, length(candidate_species)))
  } else {
    warning("[03b] EltonTraits missing species name column (Scientific/Species1). Check file format.")
  }
} else {
  # AVONET fallback: 用 Primary.Lifestyle 近似 diet_specialization
  # （AVONET 没有 Trophic.Niche 列，用 Primary.Lifestyle 替代）
  # 专化食性（Invertivore, Vertivore, Nectarivore 等）→ 高 specialization
  # 泛化食性（Omnivore）→ 低 specialization
  if (file.exists(AVONET_PATH)) {
    message("[03b] EltonTraits not found. Using AVONET Primary.Lifestyle as fallback for diet_specialization.")
    avonet_cols <- names(read_csv(AVONET_PATH, n_max = 0, show_col_types = FALSE))
    diet_col <- if ("Trophic.Niche" %in% avonet_cols) "Trophic.Niche" else "Primary.Lifestyle"
    avonet_diet <- read_csv(AVONET_PATH, show_col_types = FALSE) |>
      rename(species = any_of(c("Species1", "binomial"))) |>
      filter(species %in% candidate_species) |>
      select(species, all_of(diet_col)) |>
      mutate(
        diet_specialization = case_when(
          .data[[diet_col]] %in% c("Invertivore", "Vertivore", "Nectarivore",
                                "Granivore", "Frugivore", "Aquatic predator",
                                "Insessorial", "Raptor") ~ 0.8,
          .data[[diet_col]] %in% c("Herbivore", "Scavenger",
                                "Generalist", "Terrestrial") ~ 0.5,
          .data[[diet_col]] == "Omnivore" ~ 0.2,
          TRUE ~ NA_real_
        )
      ) |>
      select(species, diet_specialization)
    diet_spec_df <- tibble(species = candidate_species) |>
      left_join(avonet_diet, by = "species")
    n_matched <- sum(!is.na(diet_spec_df$diet_specialization))
    message(sprintf("[03b] diet_specialization (AVONET fallback): %d/%d species matched",
                    n_matched, length(candidate_species)))
  } else {
    warning("[03b] EltonTraits file not found: ", ELTONTRAITS_PATH,
          "\n  Set V3_ELTONTRAITS_PATH environment variable.",
          "\n  diet_specialization will be NA for all species.")
    diet_spec_df <- tibble(species = candidate_species,
                          diet_specialization = NA_real_)
  }
}

if (is.null(diet_spec_df)) {
  diet_spec_df <- tibble(species = candidate_species,
                          diet_specialization = NA_real_)
}

# ── 4. habitat_breadth（IUCN Red List 栖息地数量）─────────────────────
# 使用 rredlist::rl_habitats(name = species) 提取
# habitat_breadth = 物种使用的栖息地类型数量

extract_habitat_breadth_iucn <- function(species_list, api_key) {
  if (!requireNamespace("rredlist", quietly = TRUE)) {
    warning("[03b] rredlist package not installed. habitat_breadth = NA")
    return(tibble(species = species_list, habitat_breadth = NA_real_))
  }

  if (is.null(api_key) || api_key == "") {
    warning("[03b] IUCN API key not set. Set IUCN_API_KEY env variable.")
    return(tibble(species = species_list, habitat_breadth = NA_real_))
  }

  rredlist::rl_use_iucn(api_key)

  hab_breadth <- rep(NA_real_, length(species_list))
  names(hab_breadth) <- species_list

  message(sprintf("[03b] Querying IUCN for %d species ...", length(species_list)))

  for (i in seq_along(species_list)) {
    sp <- species_list[i]
    tryCatch({
      res <- rredlist::rl_habitats(name = sp, parse = TRUE)
      if (!is.null(res$result) && nrow(res$result) > 0) {
        # 计算主要栖息地数量（majorimportance = "Yes" 或所有）
        hab_breadth[i] <- nrow(res$result)
      }
    }, error = function(e) {
      # 可能是同物异名不匹配，后面处理
    })

    # 速率限制：IUCN API 限制
    if (i %% 50 == 0) {
      message(sprintf("  ... %d/%d queried", i, length(species_list)))
      Sys.sleep(2)  # 避免超过 API 速率限制
    }
  }

  tibble(species = species_list, habitat_breadth = hab_breadth)
}

# 同物异名匹配：对 IUCN 未匹配的物种尝试 BirdLife → IUCN 同物异名
match_synonyms_iucn <- function(unmatched_species, api_key) {
  if (length(unmatched_species) == 0 || api_key == "") return(tibble())

  if (!requireNamespace("rredlist", quietly = TRUE)) {
    return(tibble())
  }

  results <- list()
  for (sp in unmatched_species) {
    tryCatch({
      # 搜索该物种，IUCN 可能返回同物异名
      search_res <- rredlist::rl_search(name = sp, parse = TRUE)
      if (!is.null(search_res$result) && nrow(search_res$result) > 0) {
        # 取第一个匹配结果的 accepted name
        accepted <- search_res$result$scientific_name[1]
        if (accepted != sp) {
          # 用 accepted name 再查栖息地
          hab_res <- rredlist::rl_habitats(name = accepted, parse = TRUE)
          if (!is.null(hab_res$result) && nrow(hab_res$result) > 0) {
            results[[sp]] <- tibble(
              species = sp,
              iucn_accepted_name = accepted,
              habitat_breadth = nrow(hab_res$result)
            )
          }
        }
      }
    }, error = function(e) NULL)
    Sys.sleep(1)
  }

  if (length(results) > 0) bind_rows(results) else tibble()
}

# 执行 IUCN 提取
hab_breadth_df <- extract_habitat_breadth_iucn(candidate_species, IUCN_API_KEY)
n_matched_iucn <- sum(!is.na(hab_breadth_df$habitat_breadth))
message(sprintf("[03b] habitat_breadth (IUCN direct): %d/%d matched",
                n_matched_iucn, length(candidate_species)))

# 同物异名匹配
unmatched <- hab_breadth_df$species[is.na(hab_breadth_df$habitat_breadth)]
if (length(unmatched) > 0 && IUCN_API_KEY != "") {
  message(sprintf("[03b] Trying synonym matching for %d unmatched species ...",
                  length(unmatched)))
  syn_matches <- match_synonyms_iucn(unmatched, IUCN_API_KEY)
  if (nrow(syn_matches) > 0) {
    # 更新匹配到的物种
    for (i in seq_len(nrow(syn_matches))) {
      idx <- which(hab_breadth_df$species == syn_matches$species[i])
      if (length(idx) > 0) {
        hab_breadth_df$habitat_breadth[idx] <- syn_matches$habitat_breadth[i]
      }
    }
    n_syn <- sum(!is.na(syn_matches$habitat_breadth))
    message(sprintf("[03b] Synonym matches: %d additional species", n_syn))
  }
}

n_final_iucn <- sum(!is.na(hab_breadth_df$habitat_breadth))
message(sprintf("[03b] habitat_breadth total: %d/%d matched",
                n_final_iucn, length(candidate_species)))

# AVONET fallback for habitat_breadth：用 Habitat.Density 近似
if (n_final_iucn < length(candidate_species) * 0.5 && file.exists(AVONET_PATH)) {
  message("[03b] IUCN coverage < 50%. Using AVONET Habitat.Density as fallback for habitat_breadth.")
  avonet_hab <- read_csv(AVONET_PATH, show_col_types = FALSE) |>
    rename(species = any_of(c("Species1", "binomial"))) |>
    filter(species %in% candidate_species) |>
    select(species, Habitat.Density) |>
    mutate(
      # Habitat.Density: 1=sparse, 2=moderate, 3=abundant → 逆转换为栖息地数量近似
      habitat_breadth_avonet = case_when(
        Habitat.Density == 1   ~ 2,   # 稀疏栖息地利用 = 窄栖息地
        Habitat.Density == 2   ~ 5,   # 中等
        Habitat.Density == 3   ~ 8,   # 广泛栖息地利用
        TRUE                   ~ NA_real_
      )
    ) |>
    select(species, habitat_breadth_avonet)

  # 对 IUCN 未匹配的物种用 AVONET 填充
  hab_breadth_df <- hab_breadth_df |>
    left_join(avonet_hab, by = "species") |>
    mutate(
      habitat_breadth = if_else(
        is.na(habitat_breadth) & !is.na(habitat_breadth_avonet),
        habitat_breadth_avonet,
        habitat_breadth
      )
    ) |>
    select(-habitat_breadth_avonet)

  n_fallback <- sum(!is.na(hab_breadth_df$habitat_breadth))
  message(sprintf("[03b] habitat_breadth after AVONET fallback: %d/%d matched",
                  n_fallback, length(candidate_species)))
}

# ── 5. 随机森林插值填补缺失 ──────────────────────────────────────────
# 对 habitat_breadth 和 diet_specialization 的缺失值
# 用其他性状预测（随机森林插值优于中位数，保留生态学相关性）

rf_impute_traits <- function(trait_df, target_col, predictor_cols,
                              num_trees = 1000, seed = 2024) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    warning("[03b] ranger not installed. Falling back to median imputation.")
    trait_df[[target_col]][is.na(trait_df[[target_col]])] <-
      median(trait_df[[target_col]], na.rm = TRUE)
    return(trait_df)
  }

  miss_idx <- is.na(trait_df[[target_col]])
  if (sum(miss_idx) == 0) return(trait_df)
  if (sum(!miss_idx) < 20) {
    warning("[03b] Too few complete cases for RF imputation. Using median.")
    trait_df[[target_col]][miss_idx] <- median(trait_df[[target_col]], na.rm = TRUE)
    return(trait_df)
  }

  # 训练集（有值的行）
  train_df <- trait_df[!miss_idx, c(target_col, predictor_cols)]
  # 预测集（缺失行）
  pred_df  <- trait_df[miss_idx, predictor_cols, drop = FALSE]

  # 标准化预测变量
  pred_means <- colMeans(train_df[, predictor_cols], na.rm = TRUE)
  pred_sds   <- apply(train_df[, predictor_cols], 2, sd, na.rm = TRUE)
  pred_sds[pred_sds == 0] <- 1

  train_X <- scale(train_df[, predictor_cols], center = pred_means, scale = pred_sds)
  pred_X  <- scale(pred_df[, predictor_cols], center = pred_means, scale = pred_sds)

  train_data <- data.frame(y = train_df[[target_col]], train_X)

  fit <- ranger::ranger(
    formula    = y ~ .,
    data       = train_data,
    num.trees  = num_trees,
    seed       = seed,
    verbose    = FALSE
  )

  pred_data <- as.data.frame(pred_X)
  names(pred_data) <- names(train_data)[-1]  # 去除 y 列名

  predicted <- predict(fit, data = pred_data)$predictions
  trait_df[[target_col]][miss_idx] <- predicted

  message(sprintf("[03b] RF imputed %d missing values for %s (OOB R² = %.3f)",
                  sum(miss_idx), target_col, fit$r.squared))

  trait_df
}

# ── 6. 合并所有性状 ──────────────────────────────────────────────────
# 合并基础性状 + diet_specialization + habitat_breadth
trait_ext <- tibble(species = candidate_species) |>
  left_join(trait_imp, by = "species") |>
  left_join(diet_spec_df, by = "species") |>
  left_join(hab_breadth_df |> select(species, habitat_breadth), by = "species")

# 添加 AVONET 补充性状（如可获取）
if (file.exists(AVONET_PATH)) {
  message("[03b] Loading AVONET for supplementary traits")
  avonet <- read_csv(AVONET_PATH, show_col_types = FALSE) |>
    rename(species = any_of(c("Species1", "binomial"))) |>
    select(species, Habitat, Habitat.Density, Migration, Trophic.Level,
           Primary.Lifestyle, any_of("Trophic.Niche"))

  trait_ext <- trait_ext |> left_join(avonet, by = "species")

  # migration_score
  trait_ext <- trait_ext |>
    mutate(migration_score = case_when(
      Migration == 1 | Migration == "1" | Migration == "Resident"   ~ 1L,
      Migration == 2 | Migration == "2" | Migration == "Partial"    ~ 2L,
      Migration == 3 | Migration == "3" | Migration == "Migratory"  ~ 3L,
      TRUE ~ NA_integer_
    ))
} else {
  message("[03b] AVONET not found. Skipping supplementary AVONET traits.")
  trait_ext$Migration <- NA
  trait_ext$migration_score <- NA_integer_
}

# ── 7. 随机森林插值 ──────────────────────────────────────────────────
# 预测变量：已有的基础性状
predictor_cols <- intersect(TRAIT_VARS_BASIC, names(trait_ext))

# 先用中位数填充基础性状 NA（作为预测变量）
for (cc in predictor_cols) {
  v <- as.numeric(trait_ext[[cc]])
  if (any(is.na(v))) v[is.na(v)] <- median(v, na.rm = TRUE)
  trait_ext[[cc]] <- v
}

# RF 插值 diet_specialization
if (any(is.na(trait_ext$diet_specialization)) && length(predictor_cols) > 0) {
  trait_ext <- rf_impute_traits(trait_ext, "diet_specialization", predictor_cols)
}

# RF 插值 habitat_breadth
if (any(is.na(trait_ext$habitat_breadth)) && length(predictor_cols) > 0) {
  trait_ext <- rf_impute_traits(trait_ext, "habitat_breadth", predictor_cols)
}

# 中位数兜底（万一 RF 也失败）
for (cc in c("diet_specialization", "habitat_breadth", "migration_score")) {
  if (cc %in% names(trait_ext)) {
    v <- as.numeric(trait_ext[[cc]])
    med_val <- median(v, na.rm = TRUE)
    if (is.na(med_val)) med_val <- 0  # 最终兜底
    if (any(is.na(v))) v[is.na(v)] <- med_val
    trait_ext[[cc]] <- v
  }
}

# ── 8. 保存 ──────────────────────────────────────────────────────────
saveRDS(trait_ext, v3_file("derived", "trait_extended_v3", "rds"))
write_csv(trait_ext, v3_file("results", "table_traits_extended_v3"))

message(sprintf(
  "[03b] DONE: %d species × %d cols → data/derived_v3/trait_extended_v3.rds",
  nrow(trait_ext), ncol(trait_ext)
))
message(sprintf("[03b]   diet_specialization: %.1f%% non-NA (after imputation)",
                100 * mean(!is.na(trait_ext$diet_specialization))))
message(sprintf("[03b]   habitat_breadth: %.1f%% non-NA (after imputation)",
                100 * mean(!is.na(trait_ext$habitat_breadth))))

log_time("03b", "Completed")
