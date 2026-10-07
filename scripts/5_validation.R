# 5_validation.R: validation metrics, paired AUROC comparisons, recalibration, decision curves
# Runs after 4_models.R and uses the models as developed.

library(tidyverse)
library(pROC)
library(here)

set.seed(20260924)

load(here("data_interim", "model_cohorts.RData"))
load(here("data_interim", "sophisticated_results.RData"))
res <- sophisticated_results

y <- as.integer(res$obs_core_val)
stopifnot(identical(y, as.integer(res$obs_pre_val)), identical(pre_val$id, core_val$id))

preds <- list(
  "Primary Model"    = res$pred_l_core_val,
  "Pre-op LASSO"     = res$pred_lasso_pre_val,
  "Pre-op XGBoost"   = res$pred_xgb_pre_val,
  "Intra-op LASSO"   = res$pred_lasso_post_val,
  "Intra-op XGBoost" = res$pred_xgb_post_val
)

lp_of <- function(p) qlogis(pmin(pmax(p, 1e-4), 1 - 1e-4))

ici <- function(y, p) {
  sm <- predict(loess(y ~ p, span = 0.75, degree = 2))
  mean(abs(sm - p))
}

net_benefit <- function(y, p, t) {
  n <- length(y)
  sum(p >= t & y == 1) / n - sum(p >= t & y == 0) / n * t / (1 - t)
}

boot_ci <- function(fn, y, p, B = 1000) {
  est <- replicate(B, {
    i <- sample(length(y), replace = TRUE)
    if (length(unique(y[i])) < 2) return(NA_real_)
    fn(y[i], p[i])
  })
  unname(quantile(est, c(0.025, 0.975), na.rm = TRUE))
}

# ---- 1. VALIDATION OF THE MODELS AS DEVELOPED ----
calib_metrics <- function(y, p, model) {
  lp   <- lp_of(p)
  citl <- glm(y ~ 1, offset = lp, family = "binomial")
  slp  <- glm(y ~ lp, family = "binomial")
  roc_ <- pROC::roc(y, p, quiet = TRUE)
  auc_ci <- as.numeric(pROC::ci.auc(roc_, method = "delong"))

  o <- sum(y); e <- sum(p)
  brier_ci <- boot_ci(function(a, b) mean((a - b)^2), y, p)
  ici_ci   <- boot_ci(ici, y, p, B = 500)

  tibble(
    Model = model, n = length(y), events = o,
    AUC = auc_ci[2], AUC_LCI = auc_ci[1], AUC_UCI = auc_ci[3],
    Brier = mean((y - p)^2), Brier_LCI = brier_ci[1], Brier_UCI = brier_ci[2],
    CITL = unname(coef(citl)[1]),
    CITL_LCI = confint.default(citl)[1, 1], CITL_UCI = confint.default(citl)[1, 2],
    Slope = unname(coef(slp)[2]),
    Slope_LCI = confint.default(slp)[2, 1], Slope_UCI = confint.default(slp)[2, 2],
    OE = o / e, OE_LCI = exp(log(o / e) - 1.96 * sqrt(1 / o)), OE_UCI = exp(log(o / e) + 1.96 * sqrt(1 / o)),
    ICI = ici(y, p), ICI_LCI = ici_ci[1], ICI_UCI = ici_ci[2],
    Mean_Pred = mean(p), Max_Pred = max(p)
  )
}

val_metrics <- imap_dfr(preds, ~ calib_metrics(y, .x, .y))
print(val_metrics, width = Inf)

# ---- 2. PAIRED COMPARISON OF DISCRIMINATION AGAINST THE PRIMARY MODEL ----
roc_primary <- pROC::roc(y, preds[["Primary Model"]], quiet = TRUE)

auc_diff <- imap_dfr(preds[-1], function(p, model) {
  roc_m <- pROC::roc(y, p, quiet = TRUE)
  d <- replicate(2000, {
    i <- sample(length(y), replace = TRUE)
    if (length(unique(y[i])) < 2) return(NA_real_)
    auc_of <- function(pp) as.numeric(pROC::roc(y[i], pp[i], quiet = TRUE)$auc)
    auc_of(preds[["Primary Model"]]) - auc_of(p)
  })
  tibble(
    Model   = model,
    Diff    = as.numeric(roc_primary$auc) - as.numeric(roc_m$auc),
    LCI     = unname(quantile(d, 0.025, na.rm = TRUE)),
    UCI     = unname(quantile(d, 0.975, na.rm = TRUE)),
    p_DeLong = pROC::roc.test(roc_primary, roc_m, method = "delong", paired = TRUE)$p.value
  )
})
print(auc_diff)

# ---- 3. RECALIBRATION, CROSS-VALIDATED WITHIN THE VALIDATION COHORT ----
# 10-fold CV stratified by outcome, 20 repeats; risks are averaged over repeats
cv_recalibrate <- function(y, p, K = 10, R = 20) {
  lp  <- lp_of(p)
  out <- matrix(NA_real_, nrow = length(y), ncol = R)
  for (r in seq_len(R)) {
    folds <- integer(length(y))
    for (cls in 0:1) {
      idx <- which(y == cls)
      folds[idx] <- sample(rep_len(seq_len(K), length(idx)))
    }
    for (k in seq_len(K)) {
      tr  <- folds != k
      fit <- glm(y[tr] ~ lp[tr], family = "binomial")
      out[!tr, r] <- plogis(coef(fit)[1] + coef(fit)[2] * lp[!tr])
    }
  }
  rowMeans(out)
}

insample_recalibrate <- function(y, p) {
  predict(glm(y ~ lp_of(p), family = "binomial"), type = "response")
}

preds_recal_cv <- map(preds, ~ cv_recalibrate(y, .x))
preds_recal_in <- map(preds, ~ insample_recalibrate(y, .x))

recal_summary <- imap_dfr(preds, function(p, model) {
  p_cv <- preds_recal_cv[[model]]; p_in <- preds_recal_in[[model]]
  tibble(
    Model          = model,
    Brier_Original = mean((y - p)^2),
    Brier_InSample = mean((y - p_in)^2),
    Brier_CV       = mean((y - p_cv)^2),
    ICI_Original   = ici(y, p),
    ICI_InSample   = ici(y, p_in),
    ICI_CV         = ici(y, p_cv),
    AUC_CV         = as.numeric(pROC::roc(y, p_cv, quiet = TRUE)$auc)
  )
})
print(recal_summary, width = Inf)

# ---- 4. DECISION CURVE ANALYSIS ----
thresholds <- seq(0.01, 0.30, by = 0.005)

dca_curve <- function(pred_list, label) {
  models <- imap_dfr(pred_list, function(p, model)
    tibble(Model = model, threshold = thresholds,
           net_benefit = map_dbl(thresholds, ~ net_benefit(y, p, .x))))
  treat_all <- tibble(Model = "Treat all", threshold = thresholds,
                      net_benefit = mean(y) - (1 - mean(y)) * thresholds / (1 - thresholds))
  bind_rows(models, treat_all, tibble(Model = "Treat none", threshold = thresholds, net_benefit = 0)) |>
    mutate(Predictions = label)
}

dca_curves <- bind_rows(
  dca_curve(preds, "Original"),
  dca_curve(preds_recal_cv, "Recalibrated (cross-validated)")
)

nb_at <- c(0.05, 0.10, 0.15)
nb_table <- bind_rows(
  imap_dfr(preds, function(p, model) map_dfr(nb_at, function(t) {
    ci <- boot_ci(function(a, b) net_benefit(a, b, t), y, p)
    tibble(Model = model, Predictions = "Original", threshold = t,
           NB = net_benefit(y, p, t), NB_LCI = ci[1], NB_UCI = ci[2],
           Sensitivity = mean(p[y == 1] >= t), Specificity = mean(p[y == 0] < t),
           PPV = if (any(p >= t)) mean(y[p >= t]) else NA_real_,
           NPV = mean(1 - y[p < t]),
           Flagged_Pct = mean(p >= t))
  })),
  imap_dfr(preds_recal_cv, function(p, model) map_dfr(nb_at, function(t) {
    ci <- boot_ci(function(a, b) net_benefit(a, b, t), y, p)
    tibble(Model = model, Predictions = "Recalibrated (cross-validated)", threshold = t,
           NB = net_benefit(y, p, t), NB_LCI = ci[1], NB_UCI = ci[2],
           Sensitivity = mean(p[y == 1] >= t), Specificity = mean(p[y == 0] < t),
           PPV = if (any(p >= t)) mean(y[p >= t]) else NA_real_,
           NPV = mean(1 - y[p < t]),
           Flagged_Pct = mean(p >= t))
  })),
  tibble(Model = "Treat all", Predictions = "Original", threshold = nb_at,
         NB = mean(y) - (1 - mean(y)) * nb_at / (1 - nb_at))
)
print(nb_table |> filter(Predictions == "Original"), n = Inf)

save(val_metrics, auc_diff, recal_summary, preds_recal_cv, dca_curves, nb_table,
     file = here("data_interim", "validation_results.RData"))
message("Success: validation_results.RData saved.")
