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
set.seed(20260910)
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
    if (exists("intraop_vars")) intraop_vars else character(0)
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

CORE_VARS   <- c("age_years", "urea_t0", "cog_pre_dementia_flag")

fit_xgb <- function(X, y, nrounds = NULL, groups = NULL) {
  if (!exists("xgb_params")) stop("xgb_params not found. Run section 4 first.")
  if (is.null(nrounds)) {
    folds <- if (is.null(groups)) NULL else split(seq_along(y), grouped_folds(groups, 5))
    cv <- xgb.cv(
      params = xgb_params, data = xgb.DMatrix(X, label = y),
      nrounds = 500, nfold = 5, folds = folds, early_stopping_rounds = 20,
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

  p_hat <- predict(fit_orig, type = "response")

  boot <- future_map(seq_len(n_boot), function(i) {
    idx     <- sample(nrow(df), replace = TRUE)
    boot_df <- df[idx, ]
    if (length(unique(boot_df$delirium_7d)) < 2) return(NULL)

    fit_b  <- glm(formula_core, data = boot_df, family = "binomial")
    p_boot <- predict(fit_b, newdata = boot_df, type = "response")
    p_orig <- predict(fit_b, newdata = df,      type = "response")

    list(opt    = auc_of(boot_df$delirium_7d, p_boot) - auc_of(df$delirium_7d, p_orig),
         slope  = calib_slope(df$delirium_7d, p_orig),
         absdev = abs(p_orig - p_hat))
  }, .options = furrr_options(seed = TRUE, packages = "pROC"))

  summarise_boot(orig_auc, boot)
}

# The number of boosting rounds is re-tuned by cross-validation inside every replicate
# (fit_xgb with nrounds = NULL), so tuning is part of what the bootstrap corrects for.
get_optimism_corrected_auc_xgb <- function(df, X_orig, y_orig, fit_orig, n_boot = 1000) {
  p_hat        <- predict(fit_orig, xgb.DMatrix(X_orig))
  apparent_auc <- auc_of(y_orig, p_hat)
  message(sprintf("  -> Optimism correction xgb, nrounds re-tuned per replicate (n_boot=%d)...", n_boot))

  boot <- future_map(seq_len(n_boot), function(i) {
    idx     <- sample(seq_len(nrow(df)), replace = TRUE)
    boot_df <- df[idx, ]
    if (length(unique(boot_df$delirium_7d)) < 2) return(NULL)

    boot_obj <- make_xy(boot_df)
    X_boot   <- boot_obj$X; colnames(X_boot) <- make.names(colnames(X_boot))
    X_boot   <- align_features(X_orig, X_boot)

    fit_b  <- fit_xgb(X_boot, boot_obj$y, groups = idx)
    p_boot <- predict(fit_b, xgb.DMatrix(X_boot))
    p_orig <- predict(fit_b, xgb.DMatrix(X_orig))

    list(opt    = auc_of(boot_obj$y, p_boot) - auc_of(y_orig, p_orig),
         slope  = calib_slope(y_orig, p_orig),
         absdev = abs(p_orig - p_hat))
  }, .options = furrr_options(seed = TRUE, packages = c("xgboost", "pROC", "dplyr", "tidyr", "forcats")))

  summarise_boot(apparent_auc, boot)
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

# groups: original patient index for bootstrap samples, so that every copy of a patient
# falls in the same CV fold (otherwise duplicates leak across folds and lambda is too small)
grouped_folds <- function(groups, k) {
  ug <- unique(groups)
  sample(rep_len(seq_len(k), length(ug)))[match(groups, ug)]
}

stability_select <- function(X, y, n_reps = 100, threshold = 0.5, base_seed = 20260910, groups = NULL) {
  selected <- matrix(0L, nrow = n_reps, ncol = ncol(X), dimnames = list(NULL, colnames(X)))
  for (i in seq_len(n_reps)) {
    set.seed(base_seed + i)
    cv_i <- if (is.null(groups)) {
      cv.glmnet(X, y, family = "binomial", maxit = 1e6, thresh = 1e-7)
    } else {
      cv.glmnet(X, y, family = "binomial", foldid = grouped_folds(groups, 10), maxit = 1e6, thresh = 1e-7)
    }
    ci   <- as.matrix(coef(cv_i, s = "lambda.min"))
    nz   <- rownames(ci)[ci[, 1] != 0 & rownames(ci) != "(Intercept)"]
    selected[i, nz] <- 1L
  }
  freq <- sort(colMeans(selected), decreasing = TRUE)
  list(freq = freq, stable_vars = names(freq)[freq >= threshold])
}

refit_stable <- function(X, y, stable_vars) {
  df <- as.data.frame(X[, stable_vars, drop = FALSE])
  df$delirium_7d <- y
  glm(delirium_7d ~ ., data = df, family = "binomial")
}

calib_slope <- function(y, p) {
  lp <- qlogis(pmin(pmax(p, 1e-4), 1 - 1e-4))
  unname(coef(glm(y ~ lp, family = "binomial"))[2])
}

auc_of <- function(y, p) as.numeric(pROC::roc(y, p, quiet = TRUE)$auc)

# Harrell bootstrap summary: optimism in AUC, calibration slope of each bootstrap model
# applied to the original data (i.e. the shrinkage estimate), and prediction instability
# (mean absolute difference between bootstrap-model and final-model predictions, per patient)
summarise_boot <- function(apparent_auc, boot) {
  boot  <- Filter(Negate(is.null), boot)
  opt   <- vapply(boot, `[[`, numeric(1), "opt")
  slope <- vapply(boot, `[[`, numeric(1), "slope")
  indiv <- colMeans(do.call(rbind, lapply(boot, `[[`, "absdev")))
  tibble(Apparent_AUC           = apparent_auc,
         Optimism               = mean(opt),
         Optimism_Corrected_AUC = apparent_auc - mean(opt),
         Corrected_Slope        = mean(slope, na.rm = TRUE),
         MAPE                   = mean(indiv),
         MAPE_p95               = unname(quantile(indiv, 0.95)),
         n_boot                 = length(boot))
}

# Stability selection is repeated inside every bootstrap replicate, otherwise the
# variable-selection step escapes the optimism correction entirely.
get_optimism_corrected_auc_stab <- function(X, y, fit_orig, final_vars, n_boot = 200, n_reps = 100) {
  p_hat        <- as.numeric(predict(fit_orig, type = "response"))
  apparent_auc <- auc_of(y, p_hat)
  message(sprintf("  -> Optimism correction, stability selection repeated per replicate (n_boot=%d)...", n_boot))

  boot <- future_map(seq_len(n_boot), function(b) {
    idx <- sample(length(y), replace = TRUE)
    Xb  <- X[idx, , drop = FALSE]; yb <- y[idx]
    if (length(unique(yb)) < 2) return(NULL)

    vars_b <- stability_select(Xb, yb, n_reps, base_seed = 20260910 + 1000 * b, groups = idx)$stable_vars
    if (length(vars_b) == 0) {
      p_boot <- rep(mean(yb), length(yb))
      p_orig <- rep(mean(yb), length(y))
    } else {
      fit_b  <- refit_stable(Xb, yb, vars_b)
      p_boot <- as.numeric(predict(fit_b, type = "response"))
      p_orig <- as.numeric(predict(fit_b, newdata = as.data.frame(X[, vars_b, drop = FALSE]), type = "response"))
    }

    list(opt    = auc_of(yb, p_boot) - auc_of(y, p_orig),
         slope  = if (length(vars_b) == 0) NA_real_ else calib_slope(y, p_orig),
         absdev = abs(p_orig - p_hat),
         vars   = paste(sort(vars_b), collapse = " + "))
  }, .options = furrr_options(seed = TRUE, packages = c("glmnet", "pROC")))

  boot     <- Filter(Negate(is.null), boot)
  var_sets <- vapply(boot, `[[`, character(1), "vars")
  out      <- summarise_boot(apparent_auc, boot)
  out$Same_Set_Pct <- mean(var_sets == paste(sort(final_vars), collapse = " + "))
  attr(out, "var_sets") <- var_sets
  out
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

# 3. Reconstruct full intraoperative predictor list
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

m_x_pre  <- fit_xgb(pre_X,  pre_xy$y)
m_x_post <- fit_xgb(post_X, post_xy$y)

# LASSO variable selection is unstable with this few events, so rather than
# report one arbitrarily-seeded cv.glmnet fit, select variables by frequency across 100
# resampled fits and refit an unpenalised GLM on whatever clears the 50% threshold.
message("Running LASSO stability selection (100 resamples)...")
stab_pre  <- stability_select(pre_X,  pre_xy$y)
stab_post <- stability_select(post_X, post_xy$y)
message(sprintf("  -> Pre-op stable variables (>=50%% selection): %s", paste(stab_pre$stable_vars, collapse = ", ")))
message(sprintf("  -> Intra-op stable variables (>=50%% selection): %s", paste(stab_post$stable_vars, collapse = ", ")))

m_l_pre_stab  <- refit_stable(pre_X,  pre_xy$y,  stab_pre$stable_vars)
m_l_post_stab <- refit_stable(post_X, post_xy$y, stab_post$stable_vars)

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
pred_l_pre  <- as.numeric(predict(m_l_pre_stab,  type = "response"))
pred_x_pre  <- as.numeric(predict(m_x_pre,  xgb.DMatrix(pre_X)))
pred_l_post <- as.numeric(predict(m_l_post_stab, type = "response"))
pred_x_post <- as.numeric(predict(m_x_post, xgb.DMatrix(post_X)))
pred_l_core <- as.numeric(predict(m_l_core, type = "response"))

pred_l_pre_val  <- as.numeric(predict(m_l_pre_stab,  newdata = as.data.frame(pre_X_val[, stab_pre$stable_vars, drop = FALSE]), type = "response"))
pred_x_pre_val  <- as.numeric(predict(m_x_pre,  xgb.DMatrix(pre_X_val)))
pred_l_post_val <- as.numeric(predict(m_l_post_stab, newdata = as.data.frame(post_X_val[, stab_post$stable_vars, drop = FALSE]), type = "response"))
pred_x_post_val <- as.numeric(predict(m_x_post, xgb.DMatrix(post_X_val)))

pred_l_core_val <- as.numeric(predict(m_l_core, newdata = core_val_clean, type = "response"))
stopifnot(length(pred_l_core_val) == nrow(core_val_clean))


#  6. PERFORMANCE METRICS 
message("Performance metrics")

perf_pre_l   <- get_metrics_boot(pre_xy$y,  pred_l_pre)   |> mutate(Model = "Pre-op LASSO")
perf_pre_x   <- get_metrics_boot(pre_xy$y,  pred_x_pre)   |> mutate(Model = "Pre-op XGBoost")
perf_post_l  <- get_metrics_boot(post_xy$y, pred_l_post)  |> mutate(Model = "Intra-op LASSO")
perf_post_x  <- get_metrics_boot(post_xy$y, pred_x_post)  |> mutate(Model = "Intra-op XGBoost")
perf_core_l  <- get_metrics_boot(core_data_subset_clean$delirium_7d, pred_l_core) |> mutate(Model = "Core Model")

perf_pre_l_val  <- get_metrics_boot(pre_xy_val$y,       pred_l_pre_val)  |> mutate(Model = "Pre-op LASSO")
perf_pre_x_val  <- get_metrics_boot(pre_xy_val$y,       pred_x_pre_val)  |> mutate(Model = "Pre-op XGBoost")
perf_post_l_val <- get_metrics_boot(post_xy_val$y,      pred_l_post_val) |> mutate(Model = "Intra-op LASSO")
perf_post_x_val <- get_metrics_boot(post_xy_val$y,      pred_x_post_val) |> mutate(Model = "Intra-op XGBoost")
perf_core_l_val <- get_metrics_boot(core_val_clean$delirium_7d, pred_l_core_val) |> mutate(Model = "Core Model")

#  7. OPTIMISM CORRECTION 
message("Optimism-corrected AUC")

opt_pre_l  <- get_optimism_corrected_auc_stab(pre_X,  pre_xy$y,  m_l_pre_stab,  stab_pre$stable_vars)
opt_pre_x  <- get_optimism_corrected_auc_xgb(pre_data,  pre_X,  pre_xy$y,  m_x_pre)
opt_post_l <- get_optimism_corrected_auc_stab(post_X, post_xy$y, m_l_post_stab, stab_post$stable_vars)
opt_post_x <- get_optimism_corrected_auc_xgb(post_data, post_X, post_xy$y, m_x_post)
opt_core_l <- get_optimism_corrected_auc_glm(core_data_subset_clean)

internal_validation <- bind_rows(
  opt_core_l |> mutate(Model = "Core Model"),
  opt_pre_l  |> mutate(Model = "Pre-op LASSO"),
  opt_pre_x  |> mutate(Model = "Pre-op XGBoost"),
  opt_post_l |> mutate(Model = "Intra-op LASSO"),
  opt_post_x |> mutate(Model = "Intra-op XGBoost")
)
print(internal_validation)

# ---- 8. PERFORMANCE TABLE (LOCO-FREE VERSION) ----
message("=== Generating final performance table ===")

# 1. Compile optimism lookup for the machine learning models
opt_lookup <- tibble(
  Model = c("Pre-op LASSO", "Pre-op XGBoost", "Intra-op LASSO", "Intra-op XGBoost"),
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
  make_cal_plot(post_xy_val$y,      pred_l_post_val, "Intra-op LASSO")   +
  make_cal_plot(post_xy_val$y,      pred_x_post_val, "Intra-op XGBoost") +
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
  "LASSO (Intra-op)"   = m_l_post_roc,
  "XGBoost (Intra-op)" = m_x_post_roc
), linewidth = 1.2) +
  scale_color_manual(values = c(
    "Core Model"        = "#1b9e77",
    "LASSO (Pre-op)"    = "#d95f02",
    "XGBoost (Pre-op)"  = "#e7298a",
    "LASSO (Intra-op)"   = "#7570b3",
    "XGBoost (Intra-op)" = "#377eb8"
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
    prob_lasso_post = "Full LASSO (Intra-op)",
    prob_xgb_post   = "XGBoost (Intra-op)",
    prob_core       = "Core Model"
  )
)

p_dca_final <- plot(dca_result, smooth = TRUE) +
  theme_minimal() +
  scale_color_manual(values = c(
    "Full LASSO (Pre-op)"  = "#d95f02",
    "XGBoost (Pre-op)"     = "#e7298a",
    "Full LASSO (Intra-op)" = "#7570b3",
    "XGBoost (Intra-op)"    = "#377eb8",
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

# 1. Summarise pre-operative LASSO Model (stability-selected variables, unpenalised refit)
pre_glm_tidy  <- broom::tidy(m_l_pre_stab, conf.int = TRUE, exponentiate = FALSE)
pre_intercept <- pre_glm_tidy |> filter(term == "(Intercept)") |> pull(estimate)

coefs_pre <- pre_glm_tidy |>
  filter(term != "(Intercept)") |>
  mutate(OR = exp(estimate), LCI = exp(conf.low), UCI = exp(conf.high)) |>
  rename(Predictor = term, Beta = estimate) |>
  select(Predictor, Beta, OR, LCI, UCI) |>
  arrange(desc(abs(Beta)))

# 2. Summarise intraoperative LASSO Model (stability-selected variables, unpenalised refit)
post_glm_tidy  <- broom::tidy(m_l_post_stab, conf.int = TRUE, exponentiate = FALSE)
post_intercept <- post_glm_tidy |> filter(term == "(Intercept)") |> pull(estimate)

coefs_post <- post_glm_tidy |>
  filter(term != "(Intercept)") |>
  mutate(OR = exp(estimate), LCI = exp(conf.low), UCI = exp(conf.high)) |>
  rename(Predictor = term, Beta = estimate) |>
  select(Predictor, Beta, OR, LCI, UCI) |>
  arrange(desc(abs(Beta)))

# 3. Summarise GLM
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
  core_intercept = core_intercept, pre_intercept = pre_intercept,
  post_intercept = post_intercept, m_l_core = m_l_core,
  stab_freq_pre = stab_pre$freq, stab_freq_post = stab_post$freq,
  internal_validation = internal_validation,
  boot_var_sets_pre = attr(opt_pre_l, "var_sets"), boot_var_sets_post = attr(opt_post_l, "var_sets"),
  xgb_nrounds = c(pre = m_x_pre$niter, post = m_x_post$niter),
  core_deriv_data = core_data_subset_clean, core_val_data = core_val_clean
)

save(sophisticated_results, file = here("data_interim", "sophisticated_results.RData"))
message("Success: sophisticated_results.RData saved.")
NoSleepR::nosleep_off()


# Sample size calculation (Riley et al., BMJ 2020). A priori anticipated Nagelkerke R2 of
# 0.15 (i.e. 15% of the maximum Cox-Snell R2), not the apparent R2 of the fitted model.
library(pmsampsize)

y  <- as.integer(as.character(sophisticated_results$obs_core))
p0 <- mean(y, na.rm = TRUE)

for (k in 3:5) {
  ss <- pmsampsize(type = "b", prevalence = p0, parameters = k, nagrsquared = 0.15, shrinkage = 0.9)
  message(sprintf("  -> %d parameters: n = %d (%d events) required; available n = %d (%d events)",
                  k, ss$sample_size, ceiling(ss$sample_size * p0), length(y), sum(y)))
}
