# 4_MODELS.R  —  Model Derivation and External Validation

# ---- 1. SETUP 
library(tidyverse)
library(xgboost)
library(glmnet)
library(pROC)
library(furrr)
library(patchwork)
library(dcurves)
library(flextable)
library(broom)
library(here)

NoSleepR::nosleep_on()
n_cores <- min(4, max(1, parallel::detectCores() - 1))
plan(multisession, workers = n_cores)


# ---- 2. HELPERS 

align_features <- function(X_ref, X_new) {
  ref_cols <- colnames(X_ref)
  missing  <- setdiff(ref_cols, colnames(X_new))
  if (length(missing) > 0) {
    zero_mat <- matrix(0L, nrow = nrow(X_new), ncol = length(missing),
                       dimnames = list(NULL, missing))
    X_new <- cbind(X_new, zero_mat)
  }
  X_new[, ref_cols, drop = FALSE]
}

encode_features <- function(df) {
  df |> mutate(delirium_7d = as.integer(as.character(delirium_7d)))
}

make_xy <- function(df, outcome = "delirium_7d", medians = NULL) {
  df <- df |> mutate(delirium_7d = as.integer(as.character(delirium_7d)))
  y  <- df[[outcome]]

  feature_vars <- c(
    if (exists("t0_vars"))      t0_vars      else character(0),
    if (exists("lab_vars"))     lab_vars     else character(0),
    if (exists("intraop_vars")) intraop_vars else character(0),
    "ph_missing_flag"
  )

  df_clean <- df |>
    select(any_of(feature_vars)) |>
    select(-any_of(c(outcome, "id", "centre", "cohort_assignment"))) |>
    mutate(across(where(is.logical),   as.integer)) |>
    mutate(across(where(is.character), as.factor)) |>
    mutate(across(where(is.factor),
                  ~forcats::fct_na_value_to_level(.x, "Missing")))

  if (is.null(medians)) {
    medians <- df_clean |>
      summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)))
  }

  df_clean <- df_clean |>
    mutate(across(where(is.numeric), ~{
      fill_val <- medians[[cur_column()]]
      if (is.null(fill_val) || length(fill_val) == 0) fill_val <- median(.x, na.rm = TRUE)
      if (is.na(fill_val)) fill_val <- 0
      replace_na(.x, fill_val)
    }))

  bad_factors  <- df_clean |> select(where(is.factor))  |>
    select(where(~nlevels(droplevels(.x)) < 2)) |> names()
  bad_numerics <- df_clean |> select(where(is.numeric)) |>
    select(where(~length(unique(.x)) <= 1))     |> names()
  drop_cols <- c(bad_factors, bad_numerics)

  clinical_flag_cols <- grep("^intraop_", drop_cols, value = TRUE)
  other_drop_cols    <- setdiff(drop_cols, clinical_flag_cols)

  if (length(clinical_flag_cols) > 0)
    message(sprintf("make_xy: dropping %d uninformative clinical column(s): %s",
                    length(clinical_flag_cols), paste(clinical_flag_cols, collapse = ", ")))
  if (length(other_drop_cols) > 0)
    message(sprintf("make_xy: dropping %d other uninformative column(s): %s",
                    length(other_drop_cols), paste(other_drop_cols, collapse = ", ")))

  df_clean <- df_clean |> select(-all_of(drop_cols))

  X <- model.matrix(~ . - 1, data = df_clean)
  list(X = X, y = y, medians = medians)
}

DROPOUT_VARS <- c("ph_t0", "po2_kpa_t0", "pco2_kpa_t0", "hco3_t0", "spo2_t0")
CORE_VARS    <- c("age_years", "urea_t0", "cog_pre_dementia_flag")

fit_xgb <- function(X, y, nrounds = NULL) {
  if (!exists("xgb_params")) stop("xgb_params not found. Run section 4 first.")
  if (is.null(nrounds)) {
    cv <- xgb.cv(
      params = xgb_params, data = xgb.DMatrix(X, label = y),
      nrounds = 500, nfold = 5, early_stopping_rounds = 20,
      metrics = "logloss", verbose = 0
    )
    nrounds <- cv$best_iteration
    if (is.null(nrounds) || length(nrounds) == 0 || is.na(nrounds)) nrounds <- 200L
  }
  xgb.train(params = xgb_params, data = xgb.DMatrix(X, label = y),
            nrounds = nrounds, verbose = 0)
}

get_metrics_boot <- function(obs, pred, n_boot = 1000) {
  if (length(obs) != length(pred))
    stop(sprintf("get_metrics_boot: obs length %d != pred length %d", length(obs), length(pred)))
  keep <- !is.na(obs) & !is.na(pred)
  obs  <- obs[keep]; pred <- pred[keep]
  message(sprintf("  -> Bootstrapping 95%% CI (n=%d) on %d cores...", length(obs), n_cores))

  boot_stats <- future_map_dfr(seq_len(n_boot), function(i) {
    idx <- sample(seq_along(obs), replace = TRUE)
    o_b <- obs[idx]; p_b <- pred[idx]
    if (length(unique(o_b)) < 2) return(tibble(auc = NA_real_, brier = NA_real_))
    tibble(
      auc   = as.numeric(pROC::roc(o_b, p_b, quiet = TRUE)$auc),
      brier = mean((o_b - p_b)^2)
    )
  }, .options = furrr_options(seed = TRUE)) |> drop_na()

  tibble(
    AUC       = as.numeric(pROC::roc(obs, pred, quiet = TRUE)$auc),
    AUC_LCI   = quantile(boot_stats$auc,   0.025),
    AUC_UCI   = quantile(boot_stats$auc,   0.975),
    Brier     = mean((obs - pred)^2),
    Brier_LCI = quantile(boot_stats$brier, 0.025),
    Brier_UCI = quantile(boot_stats$brier, 0.975)
  )
}

run_loco <- function(df_full, timing_label) {
  unique_centres <- unique(df_full$centre)
  message(sprintf("  -> LOCO for %s on %d cores...", timing_label, n_cores))

  future_map_dfr(unique_centres, function(c) {
    tr_df <- df_full |> filter(centre != c)
    ts_df <- df_full |> filter(centre == c)
    if (length(unique(ts_df$delirium_7d)) < 2) return(NULL)

    tr_obj <- make_xy(tr_df)
    ts_obj <- make_xy(ts_df, medians = tr_obj$medians)
    y_test <- ts_obj$y

    X_tr <- tr_obj$X; colnames(X_tr) <- make.names(colnames(X_tr))
    X_ts <- ts_obj$X; colnames(X_ts) <- make.names(colnames(X_ts))
    X_ts <- align_features(X_tr, X_ts)

    m_l <- cv.glmnet(X_tr, tr_obj$y, family = "binomial", maxit = 1e6, thresh = 1e-7)
    p_l <- as.numeric(predict(m_l, X_ts, s = "lambda.min", type = "response")[, 1, drop = TRUE])

    m_x <- fit_xgb(X_tr, tr_obj$y)
    p_x <- predict(m_x, xgb.DMatrix(X_ts))

    tibble(
      centre    = c, timing = timing_label,
      lasso_auc = as.numeric(pROC::roc(y_test, p_l, quiet = TRUE)$auc),
      xgb_auc   = as.numeric(pROC::roc(y_test, p_x, quiet = TRUE)$auc)
    )
  }, .options = furrr_options(seed = TRUE))
}

run_loco_glm <- function(df_full) {
  df_full <- df_full |> tidyr::drop_na(all_of(CORE_VARS), delirium_7d)
  unique_centres <- unique(df_full$centre)
  message(sprintf("  -> LOCO (GLM) for Core Model on %d cores...", n_cores))

  future_map_dfr(unique_centres, function(c) {
    tr <- df_full |> filter(centre != c)
    ts <- df_full |> filter(centre == c)
    if (length(unique(ts$delirium_7d)) < 2) return(NULL)

    fit <- glm(delirium_7d ~ age_years + urea_t0 + cog_pre_dementia_flag,
               data = tr, family = "binomial")
    p   <- predict(fit, newdata = ts, type = "response")
    if (length(p) != nrow(ts)) return(NULL)

    tibble(
      centre  = c, timing = "Core Model",
      glm_auc = as.numeric(pROC::roc(ts$delirium_7d, p, quiet = TRUE)$auc),
      xgb_auc = NA_real_
    )
  }, .options = furrr_options(seed = TRUE))
}

get_optimism_corrected_auc_glm <- function(df, n_boot = 1000) {
  formula_core <- as.formula("delirium_7d ~ age_years + urea_t0 + cog_pre_dementia_flag")
  df <- df |> tidyr::drop_na(all_of(CORE_VARS), delirium_7d)

  fit_orig <- glm(formula_core, data = df, family = "binomial")
  orig_auc <- as.numeric(pROC::roc(df$delirium_7d,
                                   predict(fit_orig, type = "response"),
                                   quiet = TRUE)$auc)
  message(sprintf("  -> Optimism correction Core GLM (n=%d, n_boot=%d)...", nrow(df), n_boot))

  optimism_vals <- future_map_dbl(seq_len(n_boot), function(i) {
    idx     <- sample(nrow(df), replace = TRUE)
    boot_df <- df[idx, ]
    if (length(unique(boot_df$delirium_7d)) < 2) return(NA_real_)

    fit_b  <- glm(formula_core, data = boot_df, family = "binomial")
    p_boot <- predict(fit_b, newdata = boot_df, type = "response")
    p_orig <- predict(fit_b, newdata = df,      type = "response")

    if (length(p_boot) != nrow(boot_df) || length(p_orig) != nrow(df)) return(NA_real_)

    as.numeric(pROC::roc(boot_df$delirium_7d, p_boot, quiet = TRUE)$auc) -
      as.numeric(pROC::roc(df$delirium_7d,    p_orig, quiet = TRUE)$auc)
  }, .options = furrr_options(seed = TRUE))

  optimism <- mean(optimism_vals, na.rm = TRUE)
  tibble(Apparent_AUC           = orig_auc,
         Optimism               = optimism,
         Optimism_Corrected_AUC = orig_auc - optimism)
}

get_optimism_corrected_auc <- function(df, model_type = c("lasso", "xgb"), n_boot = NULL) {
  model_type <- match.arg(model_type)
  if (is.null(n_boot)) n_boot <- if (model_type == "xgb") 200L else 1000L

  orig_obj <- make_xy(df)
  X_orig   <- orig_obj$X; colnames(X_orig) <- make.names(colnames(X_orig))
  y_orig   <- orig_obj$y

  fixed_nrounds <- NULL

  fit_orig <- if (model_type == "lasso") {
    cv.glmnet(X_orig, y_orig, family = "binomial", maxit = 1e6, thresh = 1e-7)
  } else {
    cv <- xgb.cv(params = xgb_params, data = xgb.DMatrix(X_orig, label = y_orig),
                 nrounds = 500, nfold = 5, early_stopping_rounds = 20,
                 metrics = "logloss", verbose = 0)
    fixed_nrounds <- cv$best_iteration
    if (is.null(fixed_nrounds) || is.na(fixed_nrounds)) fixed_nrounds <- 200L
    message(sprintf("  -> XGBoost nrounds fixed at %d for bootstrap", fixed_nrounds))
    fit_xgb(X_orig, y_orig, nrounds = fixed_nrounds)
  }

  pred_orig <- if (model_type == "lasso") {
    as.numeric(predict(fit_orig, X_orig, s = "lambda.min", type = "response")[, 1, drop = TRUE])
  } else {
    predict(fit_orig, xgb.DMatrix(X_orig))
  }

  apparent_auc <- as.numeric(pROC::roc(y_orig, pred_orig, quiet = TRUE)$auc)
  message(sprintf("  -> Optimism correction %s (n_boot=%d)...", model_type, n_boot))

  optimism_vals <- future_map_dbl(seq_len(n_boot), function(i) {
    idx     <- sample(seq_len(nrow(df)), replace = TRUE)
    boot_df <- df[idx, ]
    if (length(unique(boot_df$delirium_7d)) < 2) return(NA_real_)

    boot_obj   <- make_xy(boot_df)
    X_boot_raw <- boot_obj$X; colnames(X_boot_raw) <- make.names(colnames(X_boot_raw))
    y_boot     <- boot_obj$y
    X_boot     <- align_features(X_orig, X_boot_raw)

    if (model_type == "lasso") {
      cv_b   <- glmnet::cv.glmnet(X_boot, y_boot, family = "binomial", maxit = 1e6, thresh = 1e-7)
      p_boot <- as.numeric(predict(cv_b, X_boot, s = "lambda.min", type = "response")[, 1, drop = TRUE])
      p_orig <- as.numeric(predict(cv_b, X_orig, s = "lambda.min", type = "response")[, 1, drop = TRUE])
    } else {
      fit_b  <- fit_xgb(X_boot, y_boot, nrounds = fixed_nrounds)
      p_boot <- predict(fit_b, xgb.DMatrix(X_boot))
      p_orig <- predict(fit_b, xgb.DMatrix(X_orig))
    }

    as.numeric(pROC::roc(y_boot, p_boot, quiet = TRUE)$auc) -
      as.numeric(pROC::roc(y_orig, p_orig, quiet = TRUE)$auc)
  }, .options = furrr_options(seed = TRUE))

  optimism <- mean(optimism_vals, na.rm = TRUE)
  tibble(Apparent_AUC           = apparent_auc,
         Optimism               = optimism,
         Optimism_Corrected_AUC = apparent_auc - optimism)
}

make_cal_plot <- function(obs, pred, title) {
  tibble(obs = obs, pred = pred) |>
    filter(!is.na(obs) & !is.na(pred)) |>
    ggplot(aes(x = pred, y = obs)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_smooth(method = "loess", span = 0.75,
                colour = "steelblue", fill = "steelblue", alpha = 0.2) +
    geom_rug(alpha = 0.05, sides = "b") +
    coord_cartesian(xlim = c(0, 0.8), ylim = c(0, 0.8)) +
    theme_minimal() +
    labs(title = title, x = "Predicted Probability", y = "Observed Proportion")
}

extract_lasso_coefs <- function(model, X, y, n_boot = 1000) {
  orig_coefs <- coef(model, s = "lambda.min")
  all_names  <- rownames(orig_coefs)
  
  message("Starting parallel bootstrap loop (N = ", n_boot, ")...")
  
  boot_list <- future_map(seq_len(n_boot), function(i) {
    idx  <- sample(seq_len(nrow(X)), replace = TRUE)
    cv_b <- glmnet::cv.glmnet(
      X[idx, ], y[idx], 
      family = "binomial", 
      maxit = 1e6, 
      thresh = 1e-7
    )
    
    b_coefs_raw <- coef(cv_b, s = "lambda.min")
    b_vector <- as.vector(b_coefs_raw)
    names(b_vector) <- rownames(b_coefs_raw)
    
    v <- setNames(rep(0, length(all_names)), all_names)
    v[names(b_vector)] <- b_vector
    
    return(v)
  }, .options = furrr_options(seed = TRUE))
  
  boot_df <- as.data.frame(do.call(rbind, boot_list))
  return(boot_df)
}

# 3. COHORT PREVALENCE FILTERING
message("=== Applying Prevalence Filter ===")
load(here("data_interim", "model_cohorts.RData"))

common_ids <- intersect(cohort_preop_df$id, cohort_postop_df$id)

derivation_ids <- cohort_postop_df |> 
  filter(cohort_assignment == "Derivation" & id %in% common_ids) |> 
  pull(id)

# 1. Identify swiveled drug candidate variables 
all_colnames <- colnames(cohort_postop_df)
drug_prefix_vars <- all_colnames[grep("^intraop_", all_colnames)]
intraop_drug_candidates <- setdiff(drug_prefix_vars, "intraop_hypotension")

# 2. Filter drug variables with >= 5% prevalence in derivation cohort only
valid_drug_vars <- cohort_postop_df |>
  filter(id %in% derivation_ids) |>
  summarise(across(all_of(intraop_drug_candidates), ~ mean(.x != 0, na.rm = TRUE))) |>
  pivot_longer(cols = everything(), names_to = "var", values_to = "prev") |>
  filter(prev >= 0.05) |>
  pull(var)

# 3. Reconstruct full post-operative predictor list
non_drug_intraop <- c(
  "anaes_tiva", "anaes_volatile", "anaes_regional",
  "intraop_hypotension", "hypotension_medicated",
  "estimated_blood_loss", "duration_hrs", "blood_transfusion"
)

# Merge procedural variables with the filtered drug flags
intraop_vars <- c(non_drug_intraop, valid_drug_vars)

# 4. FIT MODELS 
message("Deriving models")
pre_xy  <- make_xy(encode_features(pre_data))
pre_X   <- pre_xy$X;  colnames(pre_X)  <- make.names(colnames(pre_X))
post_xy <- make_xy(encode_features(post_data))
post_X  <- post_xy$X; colnames(post_X) <- make.names(colnames(post_X))

prevalence <- mean(pre_xy$y)

xgb_params <- list(
  objective        = "binary:logistic",
  tree_method      = "hist",
  max_depth        = 3,
  eta              = 0.01,
  subsample        = 0.6,
  colsample_bytree = 0.6,
  min_child_weight = 5,
  gamma            = 1,
  nthread          = 1,
  base_score       = prevalence
)

m_l_pre  <- cv.glmnet(pre_X,  pre_xy$y,  family = "binomial", maxit = 1e6, thresh = 1e-7)
m_x_pre  <- fit_xgb(pre_X,  pre_xy$y)
m_l_post <- cv.glmnet(post_X, post_xy$y, family = "binomial", maxit = 1e6, thresh = 1e-7)
m_x_post <- fit_xgb(post_X, post_xy$y)

core_data_subset_clean <- core_data_subset |>
  mutate(across(all_of(CORE_VARS), ~ {
    fill_val <- pre_xy$medians[[cur_column()]]
    if (is.null(fill_val) || length(fill_val) == 0 || is.na(fill_val)) fill_val <- 0
    replace_na(.x, fill_val)
  }))

m_l_core <- glm(delirium_7d ~ age_years + urea_t0 + cog_pre_dementia_flag,
                data = core_data_subset_clean, family = "binomial")

#  4a. VALIDATION MATRICES 
message("=== Preparing aligned validation matrices ===")

pre_xy_val  <- make_xy(encode_features(pre_val),  medians = pre_xy$medians)
pre_X_val   <- pre_xy_val$X;  colnames(pre_X_val)  <- make.names(colnames(pre_X_val))
pre_X_val   <- align_features(pre_X, pre_X_val)

post_xy_val <- make_xy(encode_features(post_val), medians = post_xy$medians)
post_X_val  <- post_xy_val$X; colnames(post_X_val) <- make.names(colnames(post_X_val))
post_X_val  <- align_features(post_X, post_X_val)

core_val_clean <- core_val |>
  mutate(across(all_of(CORE_VARS), ~ {
    fill_val <- pre_xy$medians[[cur_column()]]
    if (is.null(fill_val) || length(fill_val) == 0 || is.na(fill_val)) fill_val <- 0
    replace_na(.x, fill_val)
  }))

stopifnot(all(pre_xy_val$y == post_xy_val$y))
stopifnot(all(core_val_clean$id %in% pre_val$id))

#  5. PREDICTIONS 
pred_l_pre  <- as.numeric(predict(m_l_pre,  pre_X,  s = "lambda.min", type = "response")[, 1, drop = TRUE])
pred_x_pre  <- as.numeric(predict(m_x_pre,  xgb.DMatrix(pre_X)))
pred_l_post <- as.numeric(predict(m_l_post, post_X, s = "lambda.min", type = "response")[, 1, drop = TRUE])
pred_x_post <- as.numeric(predict(m_x_post, xgb.DMatrix(post_X)))
pred_l_core <- as.numeric(predict(m_l_core, type = "response"))

pred_l_pre_val  <- as.numeric(predict(m_l_pre,  pre_X_val,  s = "lambda.min", type = "response")[, 1, drop = TRUE])
pred_x_pre_val  <- as.numeric(predict(m_x_pre,  xgb.DMatrix(pre_X_val)))
pred_l_post_val <- as.numeric(predict(m_l_post, post_X_val, s = "lambda.min", type = "response")[, 1, drop = TRUE])
pred_x_post_val <- as.numeric(predict(m_x_post, xgb.DMatrix(post_X_val)))

pred_l_core_val <- as.numeric(predict(m_l_core, newdata = core_val_clean, type = "response"))
stopifnot(length(pred_l_core_val) == nrow(core_val_clean))

#  6. PERFORMANCE METRICS 
message("Performance metrics")

perf_pre_l   <- get_metrics_boot(pre_xy$y,  pred_l_pre)   |> mutate(Model = "Pre-op LASSO")
perf_pre_x   <- get_metrics_boot(pre_xy$y,  pred_x_pre)   |> mutate(Model = "Pre-op XGBoost")
perf_post_l  <- get_metrics_boot(post_xy$y, pred_l_post)  |> mutate(Model = "Post-op LASSO")
perf_post_x  <- get_metrics_boot(post_xy$y, pred_x_post)  |> mutate(Model = "Post-op XGBoost")
perf_core_l  <- get_metrics_boot(core_data_subset_clean$delirium_7d, pred_l_core) |> mutate(Model = "Core Model")

perf_pre_l_val  <- get_metrics_boot(pre_xy_val$y,       pred_l_pre_val)  |> mutate(Model = "Pre-op LASSO")
perf_pre_x_val  <- get_metrics_boot(pre_xy_val$y,       pred_x_pre_val)  |> mutate(Model = "Pre-op XGBoost")
perf_post_l_val <- get_metrics_boot(post_xy_val$y,      pred_l_post_val) |> mutate(Model = "Post-op LASSO")
perf_post_x_val <- get_metrics_boot(post_xy_val$y,      pred_x_post_val) |> mutate(Model = "Post-op XGBoost")
perf_core_l_val <- get_metrics_boot(core_val_clean$delirium_7d, pred_l_core_val) |> mutate(Model = "Core Model")

#  7. OPTIMISM CORRECTION 
message("Optimism-corrected AUC")

opt_pre_l  <- get_optimism_corrected_auc(pre_data,  "lasso")
opt_pre_x  <- get_optimism_corrected_auc(pre_data,  "xgb", n_boot = 200)
opt_post_l <- get_optimism_corrected_auc(post_data, "lasso")
opt_post_x <- get_optimism_corrected_auc(post_data, "xgb", n_boot = 200)
opt_core_l <- get_optimism_corrected_auc_glm(core_data_subset_clean)

# ---- 8. PERFORMANCE TABLE (LOCO-FREE VERSION) ----
message("=== Generating final performance table ===")

# 1. Compile optimism lookup for the machine learning models
opt_lookup <- tibble(
  Model = c("Pre-op LASSO", "Pre-op XGBoost", "Post-op LASSO", "Post-op XGBoost"),
  Optimism_Corrected_AUC = c(
    opt_pre_l$Optimism_Corrected_AUC, opt_pre_x$Optimism_Corrected_AUC,
    opt_post_l$Optimism_Corrected_AUC, opt_post_x$Optimism_Corrected_AUC
  )
)

# 2. Bind the derivation metrics together and map optimism correction
deriv_perf_summary <- bind_rows(perf_pre_l, perf_pre_x, perf_post_l, perf_post_x) |>
  left_join(opt_lookup, by = "Model") |>
  bind_rows(
    perf_core_l |> mutate(Optimism_Corrected_AUC = opt_core_l$Optimism_Corrected_AUC)
  )

# 3. Compile validation metrics
val_perf_summary <- bind_rows(
  perf_pre_l_val, perf_pre_x_val, perf_post_l_val, perf_post_x_val, perf_core_l_val
) |> 
  select(
    Model,
    AUC_val = AUC, AUC_LCI_val = AUC_LCI, AUC_UCI_val = AUC_UCI,
    Brier_val = Brier, Brier_LCI_val = Brier_LCI, Brier_UCI_val = Brier_UCI
  )

# 4. Generate final clinical performance table
final_performance_table <- deriv_perf_summary |>
  left_join(val_perf_summary, by = "Model") |>
  mutate(
    `Derivation Apparent AUC`         = sprintf("%.3f (%.3f\u2013%.3f)", AUC, AUC_LCI, AUC_UCI),
    `Derivation Corrected (Internal)` = sprintf("%.3f", Optimism_Corrected_AUC),
    `Validation External AUC`         = sprintf("%.3f (%.3f\u2013%.3f)", AUC_val, AUC_LCI_val, AUC_UCI_val),
    `Validation Brier Score`          = sprintf("%.3f (%.3f\u2013%.3f)", Brier_val, Brier_LCI_val, Brier_UCI_val)
  ) |>
  select(
    Model, 
    `Derivation Apparent AUC`, 
    `Derivation Corrected (Internal)`,
    `Validation External AUC`, 
    `Validation Brier Score`
  )

print(final_performance_table)

#  10. CALIBRATION PLOTS 
message("Generating calibration plots...")

p_cal_all <- (
  make_cal_plot(pre_xy_val$y,  pred_l_pre_val,  "Pre-op LASSO")   +
  make_cal_plot(pre_xy_val$y,  pred_x_pre_val,  "Pre-op XGBoost") +
  plot_spacer()
) / (
  make_cal_plot(post_xy_val$y,      pred_l_post_val, "Post-op LASSO")   +
  make_cal_plot(post_xy_val$y,      pred_x_post_val, "Post-op XGBoost") +
  make_cal_plot(core_val_clean$delirium_7d, pred_l_core_val, "Core Model")
) + plot_annotation(
  title    = "Calibration: Independent External Validation (LOESS)",
  subtitle = "Evaluated strictly on the Platform Trial Validation Sub-cohort"
)


# 11. ROC CURVES 
message("Generating ROC curves...")

m_l_pre_roc  <- pROC::roc(pre_xy_val$y,          pred_l_pre_val,  quiet = TRUE)
m_x_pre_roc  <- pROC::roc(pre_xy_val$y,          pred_x_pre_val,  quiet = TRUE)
m_l_post_roc <- pROC::roc(post_xy_val$y,          pred_l_post_val, quiet = TRUE)
m_x_post_roc <- pROC::roc(post_xy_val$y,          pred_x_post_val, quiet = TRUE)
m_l_core_roc <- pROC::roc(core_val_clean$delirium_7d,   pred_l_core_val, quiet = TRUE)

p_roc_internal <- pROC::ggroc(list(
  "Core Model"        = m_l_core_roc,
  "LASSO (Pre-op)"    = m_l_pre_roc,
  "XGBoost (Pre-op)"  = m_x_pre_roc,
  "LASSO (Post-op)"   = m_l_post_roc,
  "XGBoost (Post-op)" = m_x_post_roc
), linewidth = 1.2) +
  scale_color_manual(values = c(
    "Core Model"        = "#1b9e77",
    "LASSO (Pre-op)"    = "#d95f02",
    "XGBoost (Pre-op)"  = "#e7298a",
    "LASSO (Post-op)"   = "#7570b3",
    "XGBoost (Post-op)" = "#377eb8"
  )) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", alpha = 0.4) +
  theme_minimal() +
  labs(title = "Model Discrimination: Independent External Validation",
       x = "1 - Specificity", y = "Sensitivity", colour = "Model Pipeline") +
  theme(legend.position = "bottom")


# 12. DECISION CURVE ANALYSIS 
message("Running Decision Curve Analysis...")

stopifnot(length(pre_xy_val$y) == length(post_xy_val$y))
core_val_idx <- which(pre_val$id %in% core_val_clean$id)

dca_df <- tibble(
  delirium_n      = pre_xy_val$y[core_val_idx],
  prob_lasso_pre  = pred_l_pre_val[core_val_idx],
  prob_xgb_pre    = pred_x_pre_val[core_val_idx],
  prob_lasso_post = pred_l_post_val[core_val_idx],
  prob_xgb_post   = pred_x_post_val[core_val_idx],
  prob_core       = pred_l_core_val
)

stopifnot(all(dca_df$delirium_n == core_val_clean$delirium_7d))

dca_result <- dcurves::dca(
  delirium_n ~ prob_lasso_pre + prob_xgb_pre + prob_lasso_post + prob_xgb_post + prob_core,
  data       = dca_df,
  thresholds = seq(0, 0.40, by = 0.01),
  label = list(
    prob_lasso_pre  = "Full LASSO (Pre-op)",
    prob_xgb_pre    = "XGBoost (Pre-op)",
    prob_lasso_post = "Full LASSO (Post-op)",
    prob_xgb_post   = "XGBoost (Post-op)",
    prob_core       = "Core Model"
  )
)

p_dca_final <- plot(dca_result, smooth = TRUE) +
  theme_minimal() +
  scale_color_manual(values = c(
    "Full LASSO (Pre-op)"  = "#d95f02",
    "XGBoost (Pre-op)"     = "#e7298a",
    "Full LASSO (Post-op)" = "#7570b3",
    "XGBoost (Post-op)"    = "#377eb8",
    "Core Model"           = "#1b9e77",
    "All"                  = "black",
    "None"                 = "gray70"
  )) +
  coord_cartesian(ylim = c(-0.02, 0.06), xlim = c(0, 0.40)) +
  labs(title    = "Clinical Utility: Independent External Validation",
       subtitle = "Net benefit across isolated platform trial validation registry",
       x = "Threshold Probability", y = "Net Benefit", color = "Model Strategy") +
  theme(legend.position = "bottom")

print(p_dca_final)


# 13. COEFFICIENTS 
message("Extracting coefficients...")

# 1. Extract raw bootstrap data frames
raw_boot_pre  <- extract_lasso_coefs(m_l_pre,  pre_X,  pre_xy$y)
raw_boot_post <- extract_lasso_coefs(m_l_post, post_X, post_xy$y)

# 2. Summarise pre-operative LASSO Model
orig_matrix_pre <- as.matrix(coef(m_l_pre, s = "lambda.min"))
coefs_pre <- tibble(
  Predictor = rownames(orig_matrix_pre),
  Beta      = orig_matrix_pre[, 1]
) |>
  filter(Predictor != "(Intercept)") |>
  filter(Beta != 0) |> 
  mutate(
    OR  = exp(Beta),
    LCI = map_dbl(Predictor, ~ exp(quantile(raw_boot_pre[[.x]], 0.025, na.rm = TRUE))),
    UCI = map_dbl(Predictor, ~ exp(quantile(raw_boot_pre[[.x]], 0.975, na.rm = TRUE)))
  ) |>
  select(Predictor, Beta, OR, LCI, UCI) |>
  arrange(desc(abs(Beta)))

# 3. Summarise post-operative LASSO Model
orig_matrix_post <- as.matrix(coef(m_l_post, s = "lambda.min"))
coefs_post <- tibble(
  Predictor = rownames(orig_matrix_post),
  Beta      = orig_matrix_post[, 1]
) |>
  filter(Predictor != "(Intercept)") |>
  filter(Beta != 0) |> 
  mutate(
    OR  = exp(Beta),
    LCI = map_dbl(Predictor, ~ exp(quantile(raw_boot_post[[.x]], 0.025, na.rm = TRUE))),
    UCI = map_dbl(Predictor, ~ exp(quantile(raw_boot_post[[.x]], 0.975, na.rm = TRUE)))
  ) |>
  select(Predictor, Beta, OR, LCI, UCI) |>
  arrange(desc(abs(Beta)))

# 4. Summarise GLM
core_glm_tidy <- broom::tidy(m_l_core, conf.int = TRUE, exponentiate = FALSE)
core_intercept <- core_glm_tidy |> filter(term == "(Intercept)") |> pull(estimate)

coefs_core <- core_glm_tidy |>
  filter(term != "(Intercept)") |>
  mutate(OR = exp(estimate), LCI = exp(conf.low), UCI = exp(conf.high)) |>
  rename(Predictor = term, Beta = estimate) |>
  select(Predictor, Beta, OR, LCI, UCI) |>
  arrange(desc(abs(Beta)))

# ---- 14. XGBOOST IMPORTANCE ---------------------------------------------
importance_pre  <- xgb.importance(model = m_x_pre)
importance_post <- xgb.importance(model = m_x_post)


# ---- 15. SAVE -
sophisticated_results <- list(
  obs_pre = pre_xy$y,  pred_lasso_pre = pred_l_pre,  pred_xgb_pre = pred_x_pre,
  obs_post = post_xy$y, pred_lasso_post = pred_l_post, pred_xgb_post = pred_x_post,
  obs_core = core_data_subset$delirium_7d, pred_l_core = pred_l_core,

  obs_pre_val = pre_xy_val$y,         pred_lasso_pre_val  = pred_l_pre_val,
  pred_xgb_pre_val  = pred_x_pre_val, obs_post_val = post_xy_val$y,
  pred_lasso_post_val = pred_l_post_val, pred_xgb_post_val = pred_x_post_val,
  obs_core_val = core_val$delirium_7d,   pred_l_core_val   = pred_l_core_val,

  p_roc_internal = p_roc_internal, p_cal_all = p_cal_all, p_dca_final = p_dca_final,
  xgb_importance_pre = importance_pre, xgb_importance_post = importance_post,
  dca_df = dca_df, final_performance_table = final_performance_table,
  coefs_pre = coefs_pre, coefs_post = coefs_post, coefs_core = coefs_core,
  core_intercept = core_intercept, m_l_core = m_l_core
)

save(sophisticated_results, file = here("data_interim", "sophisticated_results.RData"))
message("Success: sophisticated_results.RData saved.")
NoSleepR::nosleep_off()


# Sample size calculation
library(pmsampsize)

ll_full <- as.numeric(logLik(m_l_core))

y  <- as.integer(as.character(sophisticated_results$obs_core))
p0 <- mean(y, na.rm = TRUE)
n  <- length(y)

ll_null <- sum(y * log(p0) + (1 - y) * log(1 - p0))

cox_snell <- 1 - exp((ll_null - ll_full) * (2/n))

r2_max <- 1 - exp(ll_null * (2/n))

nagelkerke <- cox_snell / r2_max

pmsampsize(
  type        = "b",
  prevalence  = p0,          
  parameters  = 3,          
  nagrsquared = nagelkerke,  
  shrinkage   = 0.9
)
