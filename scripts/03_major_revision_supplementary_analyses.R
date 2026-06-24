# ============================================
# 03. Supplementary analyses for the major revision
# Purpose:
# - Run quartile mortality summaries, sequential models, subgroup analyses,
#   missing-data summaries, and model diagnostics used in the revision package
# ============================================

locate_current_script <- function() {
  frames <- sys.frames()
  if (length(frames) > 0) {
    for (i in rev(seq_along(frames))) {
      ofile <- frames[[i]]$ofile
      if (!is.null(ofile) && nzchar(ofile)) {
        return(ofile)
      }
    }
  }

  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(sub("^--file=", "", file_arg[1]))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getSourceEditorContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) {
      return(ctx$path)
    }
  }

  NA_character_
}

find_archive_root <- function(start_dir) {
  archive_dirs <- c("01_SQL\u6570\u636e\u63d0\u53d6", "02_R\u7edf\u8ba1\u590d\u73b0", "03_\u56fe\u8868", "04_\u6587\u7ae0", "05_\u9644\u4ef6")
  repo_dirs <- c("scripts", "sql")

  for (candidate_dir in unique(c(start_dir, getwd()))) {
    current_dir <- tryCatch(
      normalizePath(candidate_dir, winslash = "/", mustWork = TRUE),
      error = function(e) NA_character_
    )
    if (length(current_dir) != 1 || is.na(current_dir) || !nzchar(current_dir)) next

    repeat {
      if (all(dir.exists(file.path(current_dir, archive_dirs))) ||
          all(dir.exists(file.path(current_dir, repo_dirs)))) {
        return(current_dir)
      }
      parent_dir <- dirname(current_dir)
      if (identical(parent_dir, current_dir)) break
      current_dir <- parent_dir
    }
  }

  stop("Unable to locate the repository/archive root automatically. Run this script from the repository root, archive root, or one of their subdirectories.")
}

this_script <- locate_current_script()
start_dir <- if (!is.na(this_script) && nzchar(this_script)) dirname(this_script) else getwd()
setwd(find_archive_root(start_dir))

rm(list = ls())

suppressPackageStartupMessages({
  library(DBI)
  library(odbc)
  library(data.table)
  library(dplyr)
  library(survival)
  library(mice)
  library(splines)
})

cache_file <- "outputs_abp_map/cache/cache_main_ge18h_Model2_SOFA_minimal.rds"
output_dir <- "outputs_abp_map/revision_major"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
m_imp <- as.integer(Sys.getenv("MIMICIV_MI_M", "40"))
maxit_imp <- as.integer(Sys.getenv("MIMICIV_MI_MAXIT", "10"))
seed_imp <- as.integer(Sys.getenv("MIMICIV_MI_SEED", "42"))

stopifnot(file.exists(cache_file))
cache <- readRDS(cache_file)

dat_lm18 <- as.data.table(cache$dat_lm18)
hour_cols <- cache$hour_cols
cut_q_main <- as.numeric(cache$cut_q_main)

has_vars <- function(x, df) x[x %in% names(df)]

extract_d1_p <- function(d1_obj) {
  if (!is.null(d1_obj$p.value)) return(as.numeric(d1_obj$p.value))
  out <- capture.output(print(d1_obj))
  idx <- grep("~~", out, fixed = TRUE)
  if (length(idx) == 0) idx <- grep("p", out, ignore.case = TRUE)
  if (length(idx) == 0) return(NA_real_)
  nums_chr <- regmatches(
    out[idx[1]],
    gregexpr("[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?", out[idx[1]], perl = TRUE)
  )[[1]]
  nums <- suppressWarnings(as.numeric(nums_chr))
  if (length(nums) >= 7) nums[7] else NA_real_
}

sanitize_covars_from_imp1 <- function(imp, covars) {
  d1 <- mice::complete(imp, 1)
  cov_use <- covars[covars %in% names(d1)]
  dropped <- c()

  for (v in cov_use) {
    x <- d1[[v]]
    if (is.character(x)) x <- factor(x)

    if (is.factor(x) && nlevels(x) < 2) {
      dropped <- c(dropped, v)
      next
    }

    if (is.numeric(x) || is.integer(x)) {
      idx <- which(!is.na(x))
      if (length(idx) == 0) {
        dropped <- c(dropped, v)
        next
      }
      ref <- x[idx[1]]
      if (isTRUE(all(is.na(x) | x == ref))) dropped <- c(dropped, v)
    }
  }

  list(covars = setdiff(cov_use, unique(dropped)), dropped = unique(dropped))
}

run_mi_hours_only <- function(df, covars, hour_cols, m = 40, maxit = 10, seed = 42, max_retry_drop_na = 2) {
  keep_cols <- unique(c("stay_id", "time_lm_days", "event_lm", covars, hour_cols))
  dat0 <- as.data.frame(df[, ..keep_cols])

  ini <- mice(dat0, maxit = 0, printFlag = FALSE)
  meth <- ini$method
  pred <- ini$predictorMatrix

  meth[] <- ""
  meth[hour_cols] <- "pmm"
  for (v in c("stay_id", "time_lm_days", "event_lm")) {
    if (v %in% names(meth)) meth[v] <- ""
  }

  pred[,] <- 0
  cov_for_pred <- intersect(covars, colnames(pred))
  for (h in 0:23) {
    target <- sprintf("h%02d", h)
    if (!target %in% rownames(pred)) next

    nb <- (h - 2):(h + 2)
    nb <- nb[nb >= 0 & nb <= 23]
    nb_cols <- sprintf("h%02d", nb)
    nb_cols <- nb_cols[nb_cols %in% colnames(pred)]

    pred[target, nb_cols] <- 1
    pred[target, target] <- 0
    if (length(cov_for_pred) > 0) pred[target, cov_for_pred] <- 1
    if ("time_lm_days" %in% colnames(pred)) pred[target, "time_lm_days"] <- 1
    if ("event_lm" %in% colnames(pred)) pred[target, "event_lm"] <- 1
  }

  set.seed(seed)
  imp <- mice(
    dat0,
    m = m,
    maxit = maxit,
    method = meth,
    predictorMatrix = pred,
    ridge = 1e-2,
    remove.collinear = TRUE,
    remove.constant = TRUE,
    printFlag = TRUE
  )

  bad_ids <- c()
  for (i in seq_len(m)) {
    di <- mice::complete(imp, i)
    na_row <- which(rowSums(is.na(di[, hour_cols, drop = FALSE])) > 0)
    if (length(na_row) > 0) bad_ids <- union(bad_ids, di$stay_id[na_row])
  }

  bad_ids <- unique(bad_ids)
  if (length(bad_ids) == 0) return(list(imp = imp, dropped_stay_id = character(0)))

  if (max_retry_drop_na <= 0) return(list(imp = imp, dropped_stay_id = bad_ids))

  df2 <- as.data.table(df)[!stay_id %in% bad_ids]
  res2 <- run_mi_hours_only(
    df = df2,
    covars = covars,
    hour_cols = hour_cols,
    m = m,
    maxit = maxit,
    seed = seed,
    max_retry_drop_na = max_retry_drop_na - 1
  )
  res2$dropped_stay_id <- unique(c(bad_ids, res2$dropped_stay_id))
  res2
}

cov_main <- unique(has_vars(c(
  "age", "gender", "weight",
  "charlson", "sofa",
  "itvtn_24h_vent_tag", "itvtn_24h_rrt_tag",
  "drug_24h_sedative_tag"
), dat_lm18))

cov_m1 <- unique(has_vars(c("age", "gender", "weight"), dat_lm18))
cov_m2 <- unique(has_vars(c("age", "gender", "weight", "charlson", "sofa"), dat_lm18))
cov_m3 <- unique(has_vars(c(
  "age", "gender", "weight", "charlson", "sofa",
  "itvtn_24h_vent_tag", "itvtn_24h_rrt_tag", "drug_24h_sedative_tag"
), dat_lm18))
cov_m4 <- unique(has_vars(c(
  "age", "gender", "weight", "charlson", "sofa",
  "itvtn_24h_vent_tag", "itvtn_24h_rrt_tag", "drug_24h_sedative_tag",
  "ne_equiv_mean_0_24h"
), dat_lm18))
cov_m5 <- unique(has_vars(c(
  "age", "gender", "weight", "charlson", "sofa",
  "itvtn_24h_vent_tag", "itvtn_24h_rrt_tag", "drug_24h_sedative_tag",
  "lab_fst_24h_lactate_first", "ne_equiv_mean_0_24h"
), dat_lm18))

mi_res <- run_mi_hours_only(dat_lm18, cov_main, hour_cols, m = m_imp, maxit = maxit_imp, seed = seed_imp)
imp <- mi_res$imp
dropped_ids <- mi_res$dropped_stay_id
fixed_extra <- unique(dat_lm18[, .(
  stay_id,
  icd_hypertension,
  ne_equiv_mean_0_24h,
  lab_fst_24h_lactate_first
)])

completed_list <- lapply(seq_len(imp$m), function(i) {
  dd <- as.data.table(mice::complete(imp, i))
  dd <- merge(dd, fixed_extra, by = "stay_id", all.x = TRUE)
  dd[, mean_map := rowMeans(.SD, na.rm = FALSE), .SDcols = hour_cols]
  dd[, map_q := cut(mean_map, breaks = cut_q_main, include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4"))]
  dd[, map_q := relevel(factor(map_q, levels = c("Q1", "Q2", "Q3", "Q4")), ref = "Q2")]
  dd[, ttr_lt65 := rowMeans(as.matrix(.SD) < 65), .SDcols = hour_cols]
  dd
})

imp1 <- copy(completed_list[[1]])

fit_from_data_list <- function(data_list, formula_obj, subset_expr = NULL) {
  lapply(data_list, function(dd) {
    ddf <- as.data.frame(copy(dd))
    if (!is.null(subset_expr)) {
      keep_idx <- eval(parse(text = subset_expr), envir = ddf, enclos = parent.frame())
      ddf <- ddf[keep_idx, , drop = FALSE]
    }
    coxph(formula_obj, data = ddf, ties = "efron")
  })
}

count_complete_case_model <- function(dd, vars) {
  use_vars <- vars[vars %in% names(dd)]
  ddf <- as.data.frame(copy(dd))
  ddf <- ddf[complete.cases(ddf[, use_vars, drop = FALSE]), , drop = FALSE]
  data.table(
    N = nrow(ddf),
    Events = sum(ddf$event_lm, na.rm = TRUE)
  )
}

pool_cox_fits <- function(fit_list) {
  coef_names <- lapply(fit_list, function(f) names(coef(f)))
  common_terms <- Reduce(intersect, coef_names)
  coef_mat <- sapply(fit_list, function(f) coef(f)[common_terms])
  if (is.null(dim(coef_mat))) coef_mat <- matrix(coef_mat, nrow = length(common_terms), dimnames = list(common_terms, NULL))
  vcov_list <- lapply(fit_list, function(f) vcov(f)[common_terms, common_terms, drop = FALSE])
  m <- length(fit_list)
  qbar <- rowMeans(coef_mat)
  ubar <- Reduce(`+`, vcov_list) / m
  bmat <- if (m > 1) stats::var(t(coef_mat)) else matrix(0, nrow = length(common_terms), ncol = length(common_terms))
  if (is.null(dim(bmat))) bmat <- matrix(bmat, nrow = length(common_terms), ncol = length(common_terms))
  tmat <- ubar + (1 + 1 / m) * bmat
  se <- sqrt(diag(tmat))
  z <- qbar / se
  p <- 2 * pnorm(-abs(z))
  list(
    qbar = qbar,
    tmat = tmat,
    se = se,
    p = p,
    terms = common_terms
  )
}

extract_map_terms <- function(pool_obj, model_label) {
  keep_terms <- c("map_qQ1", "map_qQ3", "map_qQ4")
  out <- data.table(term = keep_terms)
  out[, estimate := pool_obj$qbar[match(term, names(pool_obj$qbar))]]
  out[, se := pool_obj$se[match(term, names(pool_obj$qbar))]]
  out[, p.value := pool_obj$p[match(term, names(pool_obj$qbar))]]
  out[, Comparison := c("Q1 vs Q2", "Q3 vs Q2", "Q4 vs Q2")]
  out[, `HR (95% CI)` := sprintf(
    "%.2f (%.2f-%.2f)",
    exp(estimate),
    exp(estimate - 1.96 * se),
    exp(estimate + 1.96 * se)
  )]
  out[, Model := model_label]
  out[, .(Model, Comparison, `HR (95% CI)`, p.value)]
}

pooled_wald_p <- function(pool_obj, terms) {
  terms <- terms[terms %in% names(pool_obj$qbar)]
  if (length(terms) == 0) return(NA_real_)
  beta <- pool_obj$qbar[terms]
  vmat <- pool_obj$tmat[terms, terms, drop = FALSE]
  stat <- as.numeric(t(beta) %*% solve(vmat, beta))
  pchisq(stat, df = length(terms), lower.tail = FALSE)
}

form0 <- as.formula("Surv(time_lm_days, event_lm) ~ map_q")
form1 <- as.formula("Surv(time_lm_days, event_lm) ~ map_q + age + gender + weight")
form2 <- as.formula("Surv(time_lm_days, event_lm) ~ map_q + age + gender + weight + charlson + sofa")
form3 <- as.formula("Surv(time_lm_days, event_lm) ~ map_q + age + gender + weight + charlson + sofa + itvtn_24h_vent_tag + itvtn_24h_rrt_tag + drug_24h_sedative_tag")
form4 <- as.formula("Surv(time_lm_days, event_lm) ~ map_q + age + gender + weight + charlson + sofa + itvtn_24h_vent_tag + itvtn_24h_rrt_tag + drug_24h_sedative_tag + ne_equiv_mean_0_24h")
form5 <- as.formula("Surv(time_lm_days, event_lm) ~ map_q + age + gender + weight + charlson + sofa + itvtn_24h_vent_tag + itvtn_24h_rrt_tag + drug_24h_sedative_tag + lab_fst_24h_lactate_first + ne_equiv_mean_0_24h")

model0 <- pool_cox_fits(fit_from_data_list(completed_list, form0))
model1 <- pool_cox_fits(fit_from_data_list(completed_list, form1))
model2 <- pool_cox_fits(fit_from_data_list(completed_list, form2))
model3 <- pool_cox_fits(fit_from_data_list(completed_list, form3))
model4 <- pool_cox_fits(fit_from_data_list(completed_list, form4))
model5 <- pool_cox_fits(fit_from_data_list(completed_list, form5))

model_counts <- rbindlist(list(
  data.table(Model = "Model 0: crude", count_complete_case_model(imp1, c("time_lm_days", "event_lm", "map_q"))),
  data.table(Model = "Model 1: + age/sex/weight", count_complete_case_model(imp1, c("time_lm_days", "event_lm", "map_q", "age", "gender", "weight"))),
  data.table(Model = "Model 2: + Charlson/SOFA", count_complete_case_model(imp1, c("time_lm_days", "event_lm", "map_q", "age", "gender", "weight", "charlson", "sofa"))),
  data.table(Model = "Model 3: + MV/RRT/sedative", count_complete_case_model(imp1, c("time_lm_days", "event_lm", "map_q", "age", "gender", "weight", "charlson", "sofa", "itvtn_24h_vent_tag", "itvtn_24h_rrt_tag", "drug_24h_sedative_tag"))),
  data.table(Model = "Model 4: + NE-equivalent mean", count_complete_case_model(imp1, c("time_lm_days", "event_lm", "map_q", "age", "gender", "weight", "charlson", "sofa", "itvtn_24h_vent_tag", "itvtn_24h_rrt_tag", "drug_24h_sedative_tag", "ne_equiv_mean_0_24h"))),
  data.table(Model = "Model 5: + lactate + NE-equivalent mean", count_complete_case_model(imp1, c("time_lm_days", "event_lm", "map_q", "age", "gender", "weight", "charlson", "sofa", "itvtn_24h_vent_tag", "itvtn_24h_rrt_tag", "drug_24h_sedative_tag", "lab_fst_24h_lactate_first", "ne_equiv_mean_0_24h")))
), fill = TRUE)

seq_models <- rbindlist(list(
  extract_map_terms(model0, "Model 0: crude"),
  extract_map_terms(model1, "Model 1: + age/sex/weight"),
  extract_map_terms(model2, "Model 2: + Charlson/SOFA"),
  extract_map_terms(model3, "Model 3: + MV/RRT/sedative"),
  extract_map_terms(model4, "Model 4: + NE-equivalent mean"),
  extract_map_terms(model5, "Model 5: + lactate + NE-equivalent mean")
), fill = TRUE)

seq_models <- merge(seq_models, model_counts, by = "Model", all.x = TRUE, sort = FALSE)
setcolorder(seq_models, c("Model", "Comparison", "N", "Events", "HR (95% CI)", "p.value"))

quartile_deaths <- imp1[, .(
  N = .N,
  Deaths = sum(event_lm, na.rm = TRUE),
  Mortality_pct = round(100 * mean(event_lm, na.rm = TRUE), 1),
  TWA_MAP_median = round(median(mean_map, na.rm = TRUE), 2),
  TWA_MAP_IQR = sprintf("%.2f-%.2f", quantile(mean_map, 0.25, na.rm = TRUE), quantile(mean_map, 0.75, na.rm = TRUE))
), by = map_q][order(map_q)]

ttr_summary <- data.table(
  Metric = c(
    "TTR<65 median (IQR)",
    "Spearman correlation between TWA-MAP and TTR<65",
    "High TWA-MAP (Q4) with any TTR<65",
    "High TWA-MAP (Q4) with TTR<65 >= 0.25"
  ),
  Value = c(
    sprintf("%.3f (%.3f-%.3f)", median(imp1$ttr_lt65, na.rm = TRUE), quantile(imp1$ttr_lt65, 0.25, na.rm = TRUE), quantile(imp1$ttr_lt65, 0.75, na.rm = TRUE)),
    sprintf("%.3f", suppressWarnings(cor(imp1$mean_map, imp1$ttr_lt65, method = "spearman", use = "complete.obs"))),
    sprintf("%d/%d (%.1f%%)", sum(imp1$map_q == "Q4" & imp1$ttr_lt65 > 0, na.rm = TRUE), sum(imp1$map_q == "Q4", na.rm = TRUE), 100 * mean(imp1$ttr_lt65[imp1$map_q == "Q4"] > 0, na.rm = TRUE)),
    sprintf("%d/%d (%.1f%%)", sum(imp1$map_q == "Q4" & imp1$ttr_lt65 >= 0.25, na.rm = TRUE), sum(imp1$map_q == "Q4", na.rm = TRUE), 100 * mean(imp1$ttr_lt65[imp1$map_q == "Q4"] >= 0.25, na.rm = TRUE))
  )
)

missingness_hourly <- data.table(
  Hour = hour_cols,
  Missing_pct = round(100 * colMeans(is.na(dat_lm18[, hour_cols, with = FALSE])), 1)
)

dat_lm18[, mean_map_available := rowMeans(as.matrix(.SD), na.rm = TRUE), .SDcols = hour_cols]
missingness_summary <- data.table(
  Metric = c(
    "Patients in >=18h cohort before MI",
    "Patients dropped after post-imputation QC",
    "Per-patient missing hours, median (IQR)",
    "Observed available-hour mean MAP before MI, median (IQR)",
    "Imputed 24h TWA-MAP after MI, median (IQR)"
  ),
  Value = c(
    as.character(nrow(dat_lm18)),
    as.character(length(dropped_ids)),
    sprintf("%.0f (%.0f-%.0f)", median(rowSums(is.na(dat_lm18[, hour_cols, with = FALSE]))), quantile(rowSums(is.na(dat_lm18[, hour_cols, with = FALSE])), 0.25), quantile(rowSums(is.na(dat_lm18[, hour_cols, with = FALSE])), 0.75)),
    sprintf("%.2f (%.2f-%.2f)", median(dat_lm18$mean_map_available, na.rm = TRUE), quantile(dat_lm18$mean_map_available, 0.25, na.rm = TRUE), quantile(dat_lm18$mean_map_available, 0.75, na.rm = TRUE)),
    sprintf("%.2f (%.2f-%.2f)", median(imp1$mean_map, na.rm = TRUE), quantile(imp1$mean_map, 0.25, na.rm = TRUE), quantile(imp1$mean_map, 0.75, na.rm = TRUE))
  )
)

vaso_cut <- median(dat_lm18$ne_equiv_mean_0_24h[dat_lm18$ne_equiv_mean_0_24h > 0], na.rm = TRUE)
for (i in seq_along(completed_list)) {
  completed_list[[i]][, htn := factor(icd_hypertension, levels = c(0, 1), labels = c("No", "Yes"))]
  completed_list[[i]][, vaso_group := factor(ifelse(ne_equiv_mean_0_24h >= vaso_cut, "High", "None/Low"), levels = c("None/Low", "High"))]
}

form_htn_int <- as.formula("Surv(time_lm_days, event_lm) ~ map_q * htn + age + gender + weight + charlson + sofa + itvtn_24h_vent_tag + itvtn_24h_rrt_tag + drug_24h_sedative_tag")
form_vaso_int <- as.formula("Surv(time_lm_days, event_lm) ~ map_q * vaso_group + age + gender + weight + charlson + sofa + itvtn_24h_vent_tag + itvtn_24h_rrt_tag + drug_24h_sedative_tag")

htn_int <- pool_cox_fits(fit_from_data_list(completed_list, form_htn_int))
vaso_int <- pool_cox_fits(fit_from_data_list(completed_list, form_vaso_int))

htn_int_p <- pooled_wald_p(htn_int, c("map_qQ1:htnYes", "map_qQ3:htnYes", "map_qQ4:htnYes"))
vaso_int_p <- pooled_wald_p(vaso_int, c("map_qQ1:vaso_groupHigh", "map_qQ3:vaso_groupHigh", "map_qQ4:vaso_groupHigh"))

htn_yes <- pool_cox_fits(fit_from_data_list(completed_list, form3, subset_expr = "icd_hypertension == 1"))
htn_no <- pool_cox_fits(fit_from_data_list(completed_list, form3, subset_expr = "icd_hypertension == 0"))
vaso_high <- pool_cox_fits(fit_from_data_list(completed_list, form3, subset_expr = sprintf("ne_equiv_mean_0_24h >= %.10f", vaso_cut)))
vaso_nl <- pool_cox_fits(fit_from_data_list(completed_list, form3, subset_expr = sprintf("ne_equiv_mean_0_24h < %.10f", vaso_cut)))

make_stratum_table <- function(pool_obj, stratum_label) {
  out <- extract_map_terms(pool_obj, model_label = "")
  out[, Stratum := stratum_label]
  out[, Model := NULL]
  out
}

interaction_table <- rbindlist(list(
  make_stratum_table(htn_no, "Hypertension: no"),
  make_stratum_table(htn_yes, "Hypertension: yes"),
  make_stratum_table(vaso_nl, "Vasopressor exposure: none/low"),
  make_stratum_table(vaso_high, "Vasopressor exposure: high")
), fill = TRUE)
interaction_p <- data.table(
  Interaction = c("TWA-MAP quartile x hypertension", "TWA-MAP quartile x vasopressor exposure"),
  `P for interaction` = c(htn_int_p, vaso_int_p)
)

ph_fit <- coxph(
  Surv(time_lm_days, event_lm) ~ map_q + age + gender + weight + charlson + sofa + itvtn_24h_vent_tag + itvtn_24h_rrt_tag + drug_24h_sedative_tag,
  data = as.data.frame(imp1),
  ties = "efron"
)
ph_test <- as.data.table(cox.zph(ph_fit)$table, keep.rownames = "Term")

db <- Sys.getenv("MIMICIV_DSN", "mimic4_v31")
db_uid <- Sys.getenv("MIMICIV_DB_UID", "")
db_pwd <- Sys.getenv("MIMICIV_DB_PWD", "")
conn_args <- list(
  odbc::odbc(),
  dsn = db,
  database = Sys.getenv("MIMICIV_DATABASE", db),
  server = Sys.getenv("MIMICIV_DB_SERVER", "localhost"),
  port = as.integer(Sys.getenv("MIMICIV_DB_PORT", "5432"))
)
if (nzchar(db_uid)) conn_args$uid <- db_uid
if (nzchar(db_pwd)) conn_args$pwd <- db_pwd
con <- do.call(DBI::dbConnect, conn_args)

flow_counts <- as.data.table(dbGetQuery(con, "
WITH any_abp AS (
  SELECT COUNT(DISTINCT c.stay_id) AS any_abp
  FROM bdmcc.cohort_design_abp c
  JOIN mimiciv_icu.chartevents ce
    ON ce.stay_id = c.stay_id
  WHERE c.landmark_eligible = 1
    AND ce.charttime >= c.intime
    AND ce.charttime < c.intime + interval '24 hours'
    AND ce.itemid = 220052
    AND ce.valuenum IS NOT NULL
    AND ce.valuenum BETWEEN 20 AND 200
),
map_quality AS (
  SELECT stay_id, COUNT(*) FILTER (WHERE map_mean IS NOT NULL) AS n_hour_obs
  FROM bdmcc.map_hourly_final_abp
  GROUP BY stay_id
)
SELECT
  COUNT(*) FILTER (WHERE c.landmark_eligible = 1) AS landmark_eligible,
  (SELECT any_abp FROM any_abp) AS any_abp,
  COUNT(*) FILTER (WHERE c.landmark_eligible = 1 AND mq.n_hour_obs >= 18) AS ge18,
  COUNT(*) FILTER (WHERE c.landmark_eligible = 1 AND mq.n_hour_obs >= 20) AS ge20,
  COUNT(*) FILTER (WHERE c.landmark_eligible = 1 AND mq.n_hour_obs >= 24) AS ge24
FROM bdmcc.cohort_design_abp c
LEFT JOIN map_quality mq
  ON c.stay_id = mq.stay_id;
"))
flow_counts[, analyzed_ge18 := nrow(imp1)]
flow_counts[, analyzed_ge20 := 6268L]
flow_counts[, analyzed_ge24 := 1247L]

included_vs_excluded <- as.data.table(dbGetQuery(con, "
WITH map_quality AS (
  SELECT stay_id, COUNT(*) FILTER (WHERE map_mean IS NOT NULL) AS n_hour_obs
  FROM bdmcc.map_hourly_final_abp
  GROUP BY stay_id
),
base AS (
  SELECT
    p.stay_id,
    p.age,
    p.gender,
    p.weight,
    p.charlson,
    p.sofa,
    p.sapsii,
    p.icd_hypertension,
    p.icd_hf,
    p.icd_cad,
    p.icd_stroke,
    p.lab_fst_24h_lactate_first,
    p.itvtn_24h_vent_tag,
    p.itvtn_24h_rrt_tag,
    p.drug_24h_sedative_tag,
    COALESCE(mq.n_hour_obs, 0) AS n_hour_obs
  FROM bdmcc.bdmcc_population p
  JOIN bdmcc.cohort_design_abp c
    ON p.stay_id = c.stay_id
  LEFT JOIN map_quality mq
    ON p.stay_id = mq.stay_id
  WHERE c.landmark_eligible = 1
)
SELECT *
FROM base;
"))

included_vs_excluded[, included_main := fifelse(n_hour_obs >= 18, "Included (>=18h)", "Excluded (<18h or no ABP)")]
inc_exc_summary <- rbindlist(list(
  included_vs_excluded[, .(
    Group = unique(included_main),
    N = .N,
    Age = sprintf("%.1f (%.1f-%.1f)", median(age, na.rm = TRUE), quantile(age, 0.25, na.rm = TRUE), quantile(age, 0.75, na.rm = TRUE)),
    SOFA = sprintf("%.1f (%.1f-%.1f)", median(sofa, na.rm = TRUE), quantile(sofa, 0.25, na.rm = TRUE), quantile(sofa, 0.75, na.rm = TRUE)),
    SAPSII = sprintf("%.1f (%.1f-%.1f)", median(sapsii, na.rm = TRUE), quantile(sapsii, 0.25, na.rm = TRUE), quantile(sapsii, 0.75, na.rm = TRUE)),
    Lactate = sprintf("%.2f (%.2f-%.2f)", median(lab_fst_24h_lactate_first, na.rm = TRUE), quantile(lab_fst_24h_lactate_first, 0.25, na.rm = TRUE), quantile(lab_fst_24h_lactate_first, 0.75, na.rm = TRUE)),
    Hypertension_pct = round(100 * mean(icd_hypertension == 1, na.rm = TRUE), 1),
    Stroke_pct = round(100 * mean(icd_stroke == 1, na.rm = TRUE), 1),
    MV_pct = round(100 * mean(itvtn_24h_vent_tag == 1, na.rm = TRUE), 1),
    RRT_pct = round(100 * mean(itvtn_24h_rrt_tag == 1, na.rm = TRUE), 1)
  ), by = included_main]
), fill = TRUE)

fwrite(quartile_deaths, file.path(output_dir, "Table_R1_quartile_deaths.csv"))
fwrite(seq_models, file.path(output_dir, "Table_R2_sequential_models.csv"))
fwrite(ttr_summary, file.path(output_dir, "Table_R3_ttr_summary.csv"))
fwrite(missingness_hourly, file.path(output_dir, "Table_R4_hourly_missingness.csv"))
fwrite(missingness_summary, file.path(output_dir, "Table_R5_missingness_summary.csv"))
fwrite(interaction_table, file.path(output_dir, "Table_R6_stratified_models.csv"))
fwrite(interaction_p, file.path(output_dir, "Table_R7_interaction_tests.csv"))
fwrite(ph_test, file.path(output_dir, "Table_R8_ph_assumption.csv"))
fwrite(flow_counts, file.path(output_dir, "Table_R9_flow_counts.csv"))
fwrite(inc_exc_summary, file.path(output_dir, "Table_R10_included_vs_excluded.csv"))

try(dbDisconnect(con), silent = TRUE)

cat("\nRevision major analyses completed.\n")
cat("Output directory:", normalizePath(output_dir), "\n")
