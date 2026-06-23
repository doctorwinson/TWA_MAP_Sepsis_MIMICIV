# Observed-hour mean MAP sensitivity analysis without hourly MAP imputation
#
# This script repeats the primary adjusted Cox model after recomputing early MAP
# exposure from observed hourly invasive ABP-MAP values only.

required_pkgs <- c("data.table", "survival")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

cache_path <- file.path("outputs_abp_map", "cache", "cache_main_ge18h_Model2_SOFA_minimal.rds")
if (!file.exists(cache_path)) {
  stop("Missing cache file: ", cache_path, ". Run the main analysis script first.")
}

cache <- readRDS(cache_path)
dat <- as.data.table(cache$dat_lm18)
imp <- as.data.table(cache$imp1_hours)
hour_cols <- cache$hour_cols

dt <- dat[stay_id %in% imp$stay_id]
dt[, observed_hours_0_23 := rowSums(!is.na(.SD)), .SDcols = hour_cols]
dt[, map_observed_mean_0_23 := rowMeans(.SD, na.rm = TRUE), .SDcols = hour_cols]

hour_cols_02_23 <- hour_cols[3:24]
dt[, observed_hours_2_23 := rowSums(!is.na(.SD)), .SDcols = hour_cols_02_23]
dt[, map_observed_mean_2_23 := rowMeans(.SD, na.rm = TRUE), .SDcols = hour_cols_02_23]

format_hr <- function(hr, lo, hi) sprintf("%.2f (%.2f-%.2f)", hr, lo, hi)

fit_sensitivity <- function(data, exposure, label) {
  d <- copy(data)[is.finite(get(exposure))]
  cutpoints <- as.numeric(quantile(d[[exposure]], probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE, type = 7))
  cutpoints_unique <- cutpoints
  for (i in 2:length(cutpoints_unique)) {
    if (cutpoints_unique[i] <= cutpoints_unique[i - 1]) cutpoints_unique[i] <- cutpoints_unique[i - 1] + 1e-8
  }
  d[, map_q_sens := cut(get(exposure), breaks = cutpoints_unique, include.lowest = TRUE, right = FALSE, labels = c("Q1", "Q2", "Q3", "Q4"))]
  d[is.na(map_q_sens) & get(exposure) == max(cutpoints), map_q_sens := "Q4"]
  d[, map_q_sens := relevel(factor(map_q_sens), ref = "Q2")]

  model_vars <- c("time_lm_days", "event_lm", "map_q_sens", "age", "gender", "weight", "charlson", "sofa", "itvtn_24h_vent_tag", "itvtn_24h_rrt_tag", "drug_24h_sedative_tag")
  d_model <- d[complete.cases(d[, ..model_vars])]
  fit <- coxph(Surv(time_lm_days, event_lm) ~ map_q_sens + age + gender + weight + charlson + sofa + itvtn_24h_vent_tag + itvtn_24h_rrt_tag + drug_24h_sedative_tag, data = d_model, ties = "efron")
  sm <- summary(fit)
  quartile_deaths <- d_model[, .(Deaths = sum(event_lm == 1), N = .N), by = map_q_sens][order(map_q_sens)]

  out <- data.table(
    Analysis = label,
    Comparison = c("Q1 vs Q2", "Q3 vs Q2", "Q4 vs Q2"),
    `TWA-MAP range (mmHg)` = c(sprintf("%.2f to <%.2f", cutpoints[1], cutpoints[2]), sprintf("%.2f to <%.2f", cutpoints[3], cutpoints[4]), sprintf(">=%.2f to %.2f", cutpoints[4], cutpoints[5])),
    N = nrow(d_model),
    Events = sum(d_model$event_lm == 1),
    `Deaths by quartile` = paste0(quartile_deaths$map_q_sens, ": ", quartile_deaths$Deaths, "/", quartile_deaths$N, collapse = "; "),
    `HR (95% CI)` = c(
      format_hr(sm$conf.int["map_q_sensQ1", "exp(coef)"], sm$conf.int["map_q_sensQ1", "lower .95"], sm$conf.int["map_q_sensQ1", "upper .95"]),
      format_hr(sm$conf.int["map_q_sensQ3", "exp(coef)"], sm$conf.int["map_q_sensQ3", "lower .95"], sm$conf.int["map_q_sensQ3", "upper .95"]),
      format_hr(sm$conf.int["map_q_sensQ4", "exp(coef)"], sm$conf.int["map_q_sensQ4", "lower .95"], sm$conf.int["map_q_sensQ4", "upper .95"])
    ),
    P = c(sm$coefficients["map_q_sensQ1", "Pr(>|z|)"], sm$coefficients["map_q_sensQ3", "Pr(>|z|)"], sm$coefficients["map_q_sensQ4", "Pr(>|z|)"])
  )
  out[, P := ifelse(P < 0.001, "<0.001", sprintf("%.3f", P))]
  out
}

result <- rbindlist(list(
  fit_sensitivity(dt, "map_observed_mean_0_23", "Observed available-hour mean MAP, h00-h23, no imputation"),
  fit_sensitivity(dt[observed_hours_2_23 >= 16], "map_observed_mean_2_23", "Observed available-hour mean MAP, h02-h23, no h00-h01")
), use.names = TRUE)

dir.create(file.path("outputs_abp_map", "revision_major"), recursive = TRUE, showWarnings = FALSE)
fwrite(result, file.path("outputs_abp_map", "revision_major", "Table_R15_observed_mean_no_MI_sensitivity.csv"))
print(result)
