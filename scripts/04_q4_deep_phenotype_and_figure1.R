# ============================================
# 04. Q4 phenotype analysis and Figure 1
# Purpose:
# - Run the Q4 deep-phenotype analysis
# - Fit the phenotype-augmented model
# - Rebuild the final Figure 1 flow diagram
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
  library(ggplot2)
})

cache_file <- "outputs_abp_map/cache/cache_main_ge18h_Model2_SOFA_minimal.rds"
output_dir <- "outputs_abp_map/revision_major"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
m_imp <- as.integer(Sys.getenv("MIMICIV_MI_M", "40"))
maxit_imp <- as.integer(Sys.getenv("MIMICIV_MI_MAXIT", "10"))
seed_imp <- as.integer(Sys.getenv("MIMICIV_MI_SEED", "42"))
build_figure1 <- identical(Sys.getenv("MIMICIV_BUILD_FIGURE1", "0"), "1")

stopifnot(file.exists(cache_file))
cache <- readRDS(cache_file)

dat_lm18 <- as.data.table(cache$dat_lm18)
imp1_hours <- as.data.table(cache$imp1_hours)
hour_cols <- cache$hour_cols
cut_q_main <- as.numeric(cache$cut_q_main)

imp1_hours[, mean_map := rowMeans(.SD, na.rm = FALSE), .SDcols = hour_cols]
imp1_hours[, map_q := cut(
  mean_map,
  breaks = cut_q_main,
  include.lowest = TRUE,
  labels = c("Q1", "Q2", "Q3", "Q4")
)]
imp1_hours[, map_q := factor(map_q, levels = c("Q1", "Q2", "Q3", "Q4"))]

main_dt <- merge(dat_lm18, imp1_hours[, .(stay_id, mean_map, map_q)], by = "stay_id", all.x = TRUE)
main_dt <- main_dt[!is.na(map_q)]

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

phenotype_dt <- as.data.table(dbGetQuery(con, "
WITH main_ids AS (
  SELECT stay_id
  FROM bdmcc.base_landmark_surv_abp
)
SELECT
  p.stay_id,
  p.hadm_id,
  p.icd_hypertension,
  p.icd_stroke,
  p.icd_hf,
  p.icd_cad,
  p.itvtn_24h_cabg_tag,
  p.itvtn_24h_pci_tag,
  p.itvtn_24h_iabp_tag,
  p.itvtn_24h_picco_tag,
  p.itvtn_24h_nicom_tag,
  p.drug_24h_ne_tag,
  p.drug_24h_avp_tag,
  p.drug_24h_vaso_tag,
  i.first_careunit,
  a.admission_type
FROM main_ids m
JOIN bdmcc.bdmcc_population p
  ON m.stay_id = p.stay_id
LEFT JOIN mimiciv_icu.icustays i
  ON p.stay_id = i.stay_id
LEFT JOIN mimiciv_hosp.admissions a
  ON p.hadm_id = a.hadm_id;
"))

flow_counts <- as.data.table(dbGetQuery(con, "
WITH s0 AS (
  SELECT p.*
  FROM bdmcc.bdmcc_population p
  WHERE p.crtr_sepsis3 = 1
),
s1 AS (
  SELECT * FROM s0 WHERE icu_subject_order = 1
),
s2 AS (
  SELECT * FROM s1 WHERE (outtime - intime) >= interval '24 hours'
),
s3 AS (
  SELECT * FROM s2 WHERE age >= 18
),
s4 AS (
  SELECT * FROM s3 WHERE icd_malignancy = 0
),
s5 AS (
  SELECT * FROM s4 WHERE icd_pregnancy = 0
),
design AS (
  SELECT *
  FROM bdmcc.cohort_design_abp
),
abp_raw AS (
  SELECT DISTINCT s5.stay_id
  FROM design s5
  JOIN mimiciv_icu.chartevents ce
    ON ce.stay_id = s5.stay_id
  WHERE ce.charttime >= s5.intime
    AND ce.charttime < s5.intime + interval '24 hours'
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
  (SELECT COUNT(*) FROM s0) AS total_sepsis_icu,
  (SELECT COUNT(*) FROM s1) AS first_icu,
  (SELECT COUNT(*) FROM s2) AS los_ge24,
  (SELECT COUNT(*) FROM s3) AS adults_ge18,
  (SELECT COUNT(*) FROM s4) AS no_malignancy,
  (SELECT COUNT(*) FROM s5) AS pre_landmark_design,
  (SELECT COUNT(*) FROM design) AS landmark_eligible,
  (SELECT COUNT(*) FROM design s LEFT JOIN abp_raw a ON s.stay_id = a.stay_id WHERE a.stay_id IS NOT NULL) AS any_abp,
  (SELECT COUNT(*) FROM design s LEFT JOIN abp_raw a ON s.stay_id = a.stay_id WHERE a.stay_id IS NULL) AS no_abp,
  (SELECT COUNT(*) FROM design s JOIN abp_raw a ON s.stay_id = a.stay_id LEFT JOIN map_quality mq ON s.stay_id = mq.stay_id WHERE COALESCE(mq.n_hour_obs, 0) < 18) AS abp_lt18,
  (SELECT COUNT(*) FROM design s LEFT JOIN map_quality mq ON s.stay_id = mq.stay_id WHERE COALESCE(mq.n_hour_obs, 0) >= 18) AS ge18,
  (SELECT COUNT(*) FROM design s LEFT JOIN map_quality mq ON s.stay_id = mq.stay_id WHERE COALESCE(mq.n_hour_obs, 0) >= 20) AS ge20,
  (SELECT COUNT(*) FROM design s LEFT JOIN map_quality mq ON s.stay_id = mq.stay_id WHERE COALESCE(mq.n_hour_obs, 0) >= 24) AS ge24;
"))
flow_counts[, analyzed_ge18 := nrow(main_dt)]
flow_counts[, analyzed_ge20 := 6268L]
flow_counts[, analyzed_ge24 := 1247L]
flow_counts[, exclude_multi_icu := total_sepsis_icu - first_icu]
flow_counts[, exclude_los_lt24 := first_icu - los_ge24]
flow_counts[, exclude_malignancy := adults_ge18 - no_malignancy]
flow_counts[, exclude_pregnancy := no_malignancy - pre_landmark_design]
flow_counts[, exclude_landmark_qc := pre_landmark_design - landmark_eligible]
flow_counts[, exclude_post_mi_qc := ge18 - analyzed_ge18]

classify_careunit <- function(x) {
  dplyr::case_when(
    x %in% c(
      "Cardiac Vascular Intensive Care Unit (CVICU)",
      "Coronary Care Unit (CCU)"
    ) ~ "Cardiac ICU",
    x %in% c(
      "Medical Intensive Care Unit (MICU)",
      "Medical/Surgical Intensive Care Unit (MICU/SICU)"
    ) ~ "Medical ICU",
    x %in% c(
      "Surgical Intensive Care Unit (SICU)",
      "Trauma SICU (TSICU)",
      "PACU",
      "Surgery/Vascular/Intermediate"
    ) ~ "Surgical/Trauma ICU",
    x %in% c(
      "Neuro Intermediate",
      "Neuro Surgical Intensive Care Unit (Neuro SICU)",
      "Neuro Stepdown",
      "Neurology"
    ) ~ "Neuro ICU/Stepdown",
    TRUE ~ "Other"
  )
}

classify_admission <- function(x) {
  dplyr::case_when(
    x %in% c("EW EMER.", "DIRECT EMER.") ~ "Emergency/Direct emergent",
    x %in% c("URGENT") ~ "Urgent",
    x %in% c("ELECTIVE", "SURGICAL SAME DAY ADMISSION") ~ "Elective/Same-day surgery",
    x %in% c(
      "EU OBSERVATION",
      "OBSERVATION ADMIT",
      "DIRECT OBSERVATION",
      "AMBULATORY OBSERVATION"
    ) ~ "Observation",
    TRUE ~ "Other"
  )
}

phenotype_dt[, careunit_group := classify_careunit(first_careunit)]
phenotype_dt[, admission_group := classify_admission(admission_type)]
phenotype_dt[, cardio_proc_device := fifelse(
  coalesce(itvtn_24h_cabg_tag, 0L) == 1L |
    coalesce(itvtn_24h_pci_tag, 0L) == 1L |
    coalesce(itvtn_24h_iabp_tag, 0L) == 1L,
  1L,
  0L
)]
phenotype_dt[, advanced_hemo_monitor := fifelse(
  coalesce(itvtn_24h_picco_tag, 0L) == 1L |
    coalesce(itvtn_24h_nicom_tag, 0L) == 1L,
  1L,
  0L
)]
phenotype_dt[, any_vasopressor_tag := fifelse(coalesce(drug_24h_vaso_tag, 0L) == 1L, 1L, 0L)]

new_merge_vars <- setdiff(names(phenotype_dt), names(main_dt))
main_dt <- merge(
  main_dt,
  phenotype_dt[, c("stay_id", "hadm_id", new_merge_vars), with = FALSE],
  by = c("stay_id", "hadm_id"),
  all.x = TRUE
)

fmt_n_pct <- function(n_yes, n_total) {
  sprintf("%d/%d (%.1f%%)", n_yes, n_total, 100 * n_yes / pmax(n_total, 1))
}

safe_or <- function(exposed_q4, unexp_q4, exposed_q2, unexp_q2) {
  tab <- matrix(c(exposed_q4, unexp_q4, exposed_q2, unexp_q2), nrow = 2, byrow = TRUE)
  ft <- suppressWarnings(fisher.test(tab))
  list(
    or = unname(ft$estimate),
    lo = ft$conf.int[1],
    hi = ft$conf.int[2],
    p = ft$p.value
  )
}

binary_phenotypes <- c(
  "icd_hypertension",
  "icd_stroke",
  "icd_hf",
  "icd_cad",
  "itvtn_24h_cabg_tag",
  "itvtn_24h_pci_tag",
  "itvtn_24h_iabp_tag",
  "cardio_proc_device",
  "advanced_hemo_monitor",
  "any_vasopressor_tag",
  "drug_24h_ne_tag",
  "drug_24h_avp_tag"
)

binary_labels <- c(
  icd_hypertension = "Chronic hypertension",
  icd_stroke = "History of stroke",
  icd_hf = "Heart failure",
  icd_cad = "Coronary artery disease",
  itvtn_24h_cabg_tag = "CABG within 24 h",
  itvtn_24h_pci_tag = "PCI within 24 h",
  itvtn_24h_iabp_tag = "IABP within 24 h",
  cardio_proc_device = "Any CABG/PCI/IABP within 24 h",
  advanced_hemo_monitor = "PiCCO/NICOM within 24 h",
  any_vasopressor_tag = "Any vasopressor tag within 24 h",
  drug_24h_ne_tag = "Norepinephrine exposure tag within 24 h",
  drug_24h_avp_tag = "Vasopressin exposure tag within 24 h"
)

binary_table <- rbindlist(lapply(binary_phenotypes, function(v) {
  tmp <- main_dt[!is.na(get(v))]
  by_q <- tmp[, .(
    n_yes = sum(get(v) == 1, na.rm = TRUE),
    n_total = .N
  ), by = map_q][order(map_q)]
  q4 <- tmp[map_q == "Q4"]
  q2 <- tmp[map_q == "Q2"]
  or_obj <- safe_or(
    sum(q4[[v]] == 1, na.rm = TRUE),
    sum(q4[[v]] != 1, na.rm = TRUE),
    sum(q2[[v]] == 1, na.rm = TRUE),
    sum(q2[[v]] != 1, na.rm = TRUE)
  )
  data.table(
    Phenotype = binary_labels[[v]],
    Q1 = fmt_n_pct(by_q$n_yes[by_q$map_q == "Q1"], by_q$n_total[by_q$map_q == "Q1"]),
    Q2 = fmt_n_pct(by_q$n_yes[by_q$map_q == "Q2"], by_q$n_total[by_q$map_q == "Q2"]),
    Q3 = fmt_n_pct(by_q$n_yes[by_q$map_q == "Q3"], by_q$n_total[by_q$map_q == "Q3"]),
    Q4 = fmt_n_pct(by_q$n_yes[by_q$map_q == "Q4"], by_q$n_total[by_q$map_q == "Q4"]),
    `OR for Q4 vs Q2` = sprintf("%.2f (%.2f-%.2f)", or_obj$or, or_obj$lo, or_obj$hi),
    p.value = or_obj$p
  )
}), fill = TRUE)

categorical_profile <- function(df, group_var, label) {
  levs <- sort(unique(df[[group_var]]))
  rbindlist(lapply(levs, function(lev) {
    tmp <- copy(df)
    tmp[, present := fifelse(get(group_var) == lev, 1L, 0L)]
    by_q <- tmp[, .(
      n_yes = sum(present == 1L, na.rm = TRUE),
      n_total = .N
    ), by = map_q][order(map_q)]
    q4 <- tmp[map_q == "Q4"]
    q2 <- tmp[map_q == "Q2"]
    or_obj <- safe_or(
      sum(q4$present == 1L, na.rm = TRUE),
      sum(q4$present == 0L, na.rm = TRUE),
      sum(q2$present == 1L, na.rm = TRUE),
      sum(q2$present == 0L, na.rm = TRUE)
    )
    data.table(
      Domain = label,
      Level = lev,
      Q1 = fmt_n_pct(by_q$n_yes[by_q$map_q == "Q1"], by_q$n_total[by_q$map_q == "Q1"]),
      Q2 = fmt_n_pct(by_q$n_yes[by_q$map_q == "Q2"], by_q$n_total[by_q$map_q == "Q2"]),
      Q3 = fmt_n_pct(by_q$n_yes[by_q$map_q == "Q3"], by_q$n_total[by_q$map_q == "Q3"]),
      Q4 = fmt_n_pct(by_q$n_yes[by_q$map_q == "Q4"], by_q$n_total[by_q$map_q == "Q4"]),
      `OR for Q4 vs Q2` = sprintf("%.2f (%.2f-%.2f)", or_obj$or, or_obj$lo, or_obj$hi),
      p.value = or_obj$p
    )
  }), fill = TRUE)
}

careunit_table <- categorical_profile(main_dt[!is.na(careunit_group)], "careunit_group", "ICU careunit group")
admission_table <- categorical_profile(main_dt[!is.na(admission_group)], "admission_group", "Hospital admission type")
q4_profile_table <- rbindlist(list(binary_table, careunit_table, admission_table), fill = TRUE)

q4_dt <- main_dt[map_q == "Q4"]

fmt_median_iqr <- function(x, digits = 1) {
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
    median(x, na.rm = TRUE),
    quantile(x, 0.25, na.rm = TRUE),
    quantile(x, 0.75, na.rm = TRUE)
  )
}

q4_outcome_table <- data.table(
  Characteristic = c(
    "Age, years",
    "SOFA",
    "SAPS II",
    "Lactate, mmol/L",
    "NE-equivalent mean, ug/kg/min",
    "Chronic hypertension",
    "History of stroke",
    "Cardiac ICU",
    "Neuro ICU/Stepdown",
    "Elective/Same-day surgery admission",
    "Any CABG/PCI/IABP within 24 h"
  ),
  Survivors = c(
    fmt_median_iqr(q4_dt[event_lm == 0]$age, 1),
    fmt_median_iqr(q4_dt[event_lm == 0]$sofa, 1),
    fmt_median_iqr(q4_dt[event_lm == 0]$sapsii, 1),
    fmt_median_iqr(q4_dt[event_lm == 0]$lab_fst_24h_lactate_first, 2),
    fmt_median_iqr(q4_dt[event_lm == 0]$ne_equiv_mean_0_24h, 3),
    fmt_n_pct(sum(q4_dt[event_lm == 0]$icd_hypertension == 1, na.rm = TRUE), nrow(q4_dt[event_lm == 0])),
    fmt_n_pct(sum(q4_dt[event_lm == 0]$icd_stroke == 1, na.rm = TRUE), nrow(q4_dt[event_lm == 0])),
    fmt_n_pct(sum(q4_dt[event_lm == 0]$careunit_group == "Cardiac ICU", na.rm = TRUE), nrow(q4_dt[event_lm == 0])),
    fmt_n_pct(sum(q4_dt[event_lm == 0]$careunit_group == "Neuro ICU/Stepdown", na.rm = TRUE), nrow(q4_dt[event_lm == 0])),
    fmt_n_pct(sum(q4_dt[event_lm == 0]$admission_group == "Elective/Same-day surgery", na.rm = TRUE), nrow(q4_dt[event_lm == 0])),
    fmt_n_pct(sum(q4_dt[event_lm == 0]$cardio_proc_device == 1, na.rm = TRUE), nrow(q4_dt[event_lm == 0]))
  ),
  Non_survivors = c(
    fmt_median_iqr(q4_dt[event_lm == 1]$age, 1),
    fmt_median_iqr(q4_dt[event_lm == 1]$sofa, 1),
    fmt_median_iqr(q4_dt[event_lm == 1]$sapsii, 1),
    fmt_median_iqr(q4_dt[event_lm == 1]$lab_fst_24h_lactate_first, 2),
    fmt_median_iqr(q4_dt[event_lm == 1]$ne_equiv_mean_0_24h, 3),
    fmt_n_pct(sum(q4_dt[event_lm == 1]$icd_hypertension == 1, na.rm = TRUE), nrow(q4_dt[event_lm == 1])),
    fmt_n_pct(sum(q4_dt[event_lm == 1]$icd_stroke == 1, na.rm = TRUE), nrow(q4_dt[event_lm == 1])),
    fmt_n_pct(sum(q4_dt[event_lm == 1]$careunit_group == "Cardiac ICU", na.rm = TRUE), nrow(q4_dt[event_lm == 1])),
    fmt_n_pct(sum(q4_dt[event_lm == 1]$careunit_group == "Neuro ICU/Stepdown", na.rm = TRUE), nrow(q4_dt[event_lm == 1])),
    fmt_n_pct(sum(q4_dt[event_lm == 1]$admission_group == "Elective/Same-day surgery", na.rm = TRUE), nrow(q4_dt[event_lm == 1])),
    fmt_n_pct(sum(q4_dt[event_lm == 1]$cardio_proc_device == 1, na.rm = TRUE), nrow(q4_dt[event_lm == 1]))
  )
)

has_vars <- function(x, df) x[x %in% names(df)]

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

fit_from_data_list <- function(data_list, formula_obj) {
  lapply(data_list, function(dd) {
    coxph(formula_obj, data = as.data.frame(copy(dd)), ties = "efron")
  })
}

pool_cox_fits <- function(fit_list) {
  coef_names <- lapply(fit_list, function(f) names(coef(f)))
  common_terms <- Reduce(intersect, coef_names)
  coef_mat <- sapply(fit_list, function(f) coef(f)[common_terms])
  if (is.null(dim(coef_mat))) {
    coef_mat <- matrix(coef_mat, nrow = length(common_terms), dimnames = list(common_terms, NULL))
  }
  vcov_list <- lapply(fit_list, function(f) vcov(f)[common_terms, common_terms, drop = FALSE])
  m <- length(fit_list)
  qbar <- rowMeans(coef_mat)
  ubar <- Reduce(`+`, vcov_list) / m
  bmat <- if (m > 1) stats::var(t(coef_mat)) else matrix(0, nrow = length(common_terms), ncol = length(common_terms))
  if (is.null(dim(bmat))) {
    bmat <- matrix(bmat, nrow = length(common_terms), ncol = length(common_terms))
  }
  tmat <- ubar + (1 + 1 / m) * bmat
  se <- sqrt(diag(tmat))
  z <- qbar / se
  p <- 2 * pnorm(-abs(z))
  list(qbar = qbar, tmat = tmat, se = se, p = p)
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

cov_main <- unique(has_vars(c(
  "age", "gender", "weight",
  "charlson", "sofa",
  "itvtn_24h_vent_tag", "itvtn_24h_rrt_tag",
  "drug_24h_sedative_tag"
), dat_lm18))

mi_res <- run_mi_hours_only(dat_lm18, cov_main, hour_cols, m = m_imp, maxit = maxit_imp, seed = seed_imp)
imp <- mi_res$imp

completed_list <- lapply(seq_len(imp$m), function(i) {
  dd <- as.data.table(mice::complete(imp, i))
  dd[, mean_map := rowMeans(.SD, na.rm = FALSE), .SDcols = hour_cols]
  dd[, map_q := cut(mean_map, breaks = cut_q_main, include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4"))]
  dd[, map_q := relevel(factor(map_q, levels = c("Q1", "Q2", "Q3", "Q4")), ref = "Q2")]
  dup_vars <- intersect(c("icd_hypertension", "icd_stroke", "careunit_group", "admission_group", "cardio_proc_device"), names(dd))
  if (length(dup_vars) > 0) dd[, (dup_vars) := NULL]
  dd <- merge(dd, phenotype_dt[, .(
    stay_id,
    icd_hypertension,
    icd_stroke,
    careunit_group,
    admission_group,
    cardio_proc_device
  )], by = "stay_id", all.x = TRUE)
  dd[, icd_hypertension := factor(icd_hypertension, levels = c(0, 1), labels = c("No", "Yes"))]
  dd[, icd_stroke := factor(icd_stroke, levels = c(0, 1), labels = c("No", "Yes"))]
  dd[, cardio_proc_device := factor(cardio_proc_device, levels = c(0, 1), labels = c("No", "Yes"))]
  dd[, careunit_group := factor(
    careunit_group,
    levels = c("Medical ICU", "Cardiac ICU", "Surgical/Trauma ICU", "Neuro ICU/Stepdown", "Other")
  )]
  dd[, admission_group := factor(
    admission_group,
    levels = c("Emergency/Direct emergent", "Urgent", "Elective/Same-day surgery", "Observation", "Other")
  )]
  dd
})

form_pheno <- as.formula(paste(
  "Surv(time_lm_days, event_lm) ~ map_q + age + gender + weight + charlson + sofa +",
  "itvtn_24h_vent_tag + itvtn_24h_rrt_tag + drug_24h_sedative_tag +",
  "icd_hypertension + icd_stroke + careunit_group + admission_group + cardio_proc_device"
))

pheno_model <- pool_cox_fits(fit_from_data_list(completed_list, form_pheno))
pheno_model_table <- extract_map_terms(pheno_model, "Phenotype-augmented sensitivity model")
pheno_model_overall <- data.table(
  Metric = c("Overall P for map_q terms"),
  Value = pooled_wald_p(pheno_model, c("map_qQ1", "map_qQ3", "map_qQ4"))
)

if (build_figure1) {
flow_box <- function(xmin, xmax, ymin, ymax, label, fill, size = 0.38, text_size = 3.95, lineheight = 1.00) {
  list(
    annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
             fill = fill, color = "#24333b", linewidth = size),
    annotate("text", x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
             label = label, family = "sans", size = text_size, lineheight = lineheight, color = "#0f1d24")
  )
}

fc <- flow_counts[1]

main_nodes <- list(
  list(y = 13.6, label = sprintf("Sepsis-3 ICU stays in MIMIC-IV v3.1\nn = %s", format(fc$total_sepsis_icu, big.mark = ",")), fill = "#d7efe9"),
  list(y = 11.75, label = sprintf("First ICU admission\nn = %s", format(fc$first_icu, big.mark = ",")), fill = "#e8f5f2"),
  list(y = 9.9, label = sprintf("Adult first-ICU stays with ICU LOS \u226524 h\nn = %s", format(fc$adults_ge18, big.mark = ",")), fill = "#edf4ef"),
  list(y = 8.05, label = sprintf("No malignancy\nn = %s", format(fc$no_malignancy, big.mark = ",")), fill = "#edf4ef"),
  list(y = 6.2, label = sprintf("Pre-landmark design cohort\nn = %s", format(fc$pre_landmark_design, big.mark = ",")), fill = "#f4f0e8"),
  list(y = 4.35, label = sprintf("24-h landmark-eligible cohort\nn = %s", format(fc$landmark_eligible, big.mark = ",")), fill = "#f2eadf"),
  list(y = 2.5, label = sprintf("Any invasive ABP-MAP in first 24 h\n(itemid 220052)\nn = %s", format(fc$any_abp, big.mark = ",")), fill = "#f6ebd9"),
  list(y = 0.75, label = sprintf("\u226518 observed hourly ABP-MAP values\nn = %s", format(fc$ge18, big.mark = ",")), fill = "#f8dfc5"),
  list(y = -1.35, label = sprintf("Primary analyzed cohort\n\u226518 h + post-imputation QC pass\nn = %s", format(fc$analyzed_ge18, big.mark = ",")), fill = "#f3c7a6")
)

exclude_nodes <- list(
  list(y = 11.75, label = sprintf("Excluded repeat ICU admissions\nn = %s", format(fc$exclude_multi_icu, big.mark = ","))),
  list(y = 9.9, label = sprintf("Excluded ICU LOS <24 h\nn = %s", format(fc$exclude_los_lt24, big.mark = ","))),
  list(y = 8.05, label = sprintf("Excluded malignancy\nn = %s", format(fc$exclude_malignancy, big.mark = ","))),
  list(y = 6.2, label = sprintf("Excluded pregnancy\nn = %s", format(fc$exclude_pregnancy, big.mark = ","))),
  list(y = 4.35, label = sprintf("Not 24-h landmark-eligible / design QC\nn = %s", format(fc$exclude_landmark_qc, big.mark = ","))),
  list(y = 2.5, label = sprintf("No valid invasive ABP-MAP in first 24 h\nn = %s", format(fc$no_abp, big.mark = ","))),
  list(y = 0.75, label = sprintf("<18 observed hourly ABP-MAP values\nn = %s", format(fc$abp_lt18, big.mark = ","))),
  list(y = -1.35, label = sprintf("Post-imputation QC failure\nn = %s", format(fc$exclude_post_mi_qc, big.mark = ",")))
)

branch_nodes <- list(
  list(x = -3.35, y = -4.20, label = sprintf("Sensitivity cohort\n\u226520 observed hours\nn = %s", format(fc$analyzed_ge20, big.mark = ",")), fill = "#f0dcc9"),
  list(x = 3.35, y = -4.20, label = sprintf("Sensitivity cohort\n24/24 observed hours\nn = %s", format(fc$analyzed_ge24, big.mark = ",")), fill = "#ead3bc")
)

p <- ggplot() +
  coord_cartesian(xlim = c(-8.3, 8.4), ylim = c(-6.3, 14.4), clip = "off") +
  theme_void(base_family = "sans") +
  theme(
    plot.background = element_rect(fill = "#fbfaf6", color = NA),
    panel.background = element_rect(fill = "#fbfaf6", color = NA),
    plot.margin = margin(10, 34, 20, 34)
  )

for (i in seq_along(main_nodes)) {
  nd <- main_nodes[[i]]
  p <- p + flow_box(-2.60, 2.60, nd$y - 0.64, nd$y + 0.64, nd$label, nd$fill, text_size = 3.95, lineheight = 1.00)
}

for (i in seq_along(exclude_nodes)) {
  nd <- exclude_nodes[[i]]
  p <- p + flow_box(4.10, 8.00, nd$y - 0.55, nd$y + 0.55, nd$label, "#ffffff", size = 0.32, text_size = 3.70, lineheight = 0.99)
}

for (i in seq_along(branch_nodes)) {
  nd <- branch_nodes[[i]]
  p <- p + flow_box(nd$x - 2.20, nd$x + 2.20, nd$y - 0.64, nd$y + 0.64, nd$label, nd$fill, size = 0.32, text_size = 3.72, lineheight = 1.00)
}

for (y in c(13.6, 11.75, 9.9, 8.05, 6.2, 4.35, 2.5, 0.75)) {
  p <- p + annotate("segment", x = 0, xend = 0, y = y - 0.64, yend = y - 1.21,
                    linewidth = 0.42, color = "#3c4a52",
                    arrow = arrow(length = unit(0.14, "inches"), type = "closed"))
}

for (y in c(11.75, 9.9, 8.05, 6.2, 4.35, 2.5, 0.75, -1.35)) {
  p <- p + annotate("segment", x = 2.60, xend = 4.10, y = y, yend = y,
                    linewidth = 0.34, color = "#58656d",
                    arrow = arrow(length = unit(0.12, "inches"), type = "closed"))
}

p <- p +
  annotate("segment", x = 0, xend = -3.35, y = -2.02, yend = -3.56,
           linewidth = 0.34, color = "#55636b",
           arrow = arrow(length = unit(0.12, "inches"), type = "closed")) +
  annotate("segment", x = 0, xend = 3.35, y = -2.02, yend = -3.56,
           linewidth = 0.34, color = "#55636b",
           arrow = arrow(length = unit(0.12, "inches"), type = "closed")) +
  annotate("text", x = 0, y = -5.85,
           label = "Primary cohort required \u226518 observed hourly invasive ABP-MAP values in the first 24 h; sensitivity cohorts used \u226520 h and complete 24/24 h observation thresholds.",
           family = "sans", size = 3.45, color = "#33424a")

flow_box <- function(xmin, xmax, ymin, ymax, label, fill, size = 0.38, text_size = 3.80, lineheight = 0.98) {
  list(
    annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
             fill = fill, color = "#24333b", linewidth = size),
    annotate("text", x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
             label = label, family = "sans", size = text_size, lineheight = lineheight, color = "#0f1d24")
  )
}

ge <- "\u2265"

main_nodes <- list(
  list(y = 13.6, label = sprintf("Sepsis-3 ICU stays in MIMIC-IV v3.1\nn = %s", format(fc$total_sepsis_icu, big.mark = ",")), fill = "#d7efe9"),
  list(y = 11.75, label = sprintf("First ICU admission\nn = %s", format(fc$first_icu, big.mark = ",")), fill = "#e8f5f2"),
  list(y = 9.9, label = sprintf("Adult first-ICU stays\nwith ICU LOS %s24 h\nn = %s", ge, format(fc$adults_ge18, big.mark = ",")), fill = "#edf4ef"),
  list(y = 8.05, label = sprintf("No malignancy\nn = %s", format(fc$no_malignancy, big.mark = ",")), fill = "#edf4ef"),
  list(y = 6.2, label = sprintf("Pre-landmark design cohort\nn = %s", format(fc$pre_landmark_design, big.mark = ",")), fill = "#f4f0e8"),
  list(y = 4.35, label = sprintf("24-h landmark-eligible cohort\nn = %s", format(fc$landmark_eligible, big.mark = ",")), fill = "#f2eadf"),
  list(y = 2.5, label = sprintf("Any invasive ABP-MAP\nin first 24 h\n(itemid 220052)\nn = %s", format(fc$any_abp, big.mark = ",")), fill = "#f6ebd9"),
  list(y = 0.75, label = sprintf("%s18/24 observed\nhourly ABP-MAP\nvalues\nn = %s", ge, format(fc$ge18, big.mark = ",")), fill = "#f8dfc5"),
  list(y = -1.50, label = sprintf("Primary analyzed cohort\n%s18/24 h observed +\npost-imputation QC pass\nn = %s", ge, format(fc$analyzed_ge18, big.mark = ",")), fill = "#f3c7a6")
)

exclude_nodes <- list(
  list(y = 11.75, label = sprintf("Excluded repeat ICU admissions\nn = %s", format(fc$exclude_multi_icu, big.mark = ","))),
  list(y = 9.9, label = sprintf("Excluded ICU LOS <24 h\nn = %s", format(fc$exclude_los_lt24, big.mark = ","))),
  list(y = 8.05, label = sprintf("Excluded malignancy\nn = %s", format(fc$exclude_malignancy, big.mark = ","))),
  list(y = 6.2, label = sprintf("Excluded pregnancy\nn = %s", format(fc$exclude_pregnancy, big.mark = ","))),
  list(y = 4.35, label = sprintf("Excluded: not 24-h\nlandmark-eligible or\ndesign QC failure\nn = %s", format(fc$exclude_landmark_qc, big.mark = ","))),
  list(y = 2.5, label = sprintf("No valid invasive\nABP-MAP in first\n24 h\nn = %s", format(fc$no_abp, big.mark = ","))),
  list(y = 0.75, label = sprintf("<18/24 observed\nhourly ABP-MAP\nvalues\nn = %s", format(fc$abp_lt18, big.mark = ","))),
  list(y = -1.50, label = sprintf("Post-imputation QC failure\nn = %s", format(fc$exclude_post_mi_qc, big.mark = ",")))
)

branch_nodes <- list(
  list(x = -3.45, y = -4.45, label = sprintf("Sensitivity cohort\n%s20/24 observed hours\nn = %s", ge, format(fc$analyzed_ge20, big.mark = ",")), fill = "#f0dcc9"),
  list(x = 3.45, y = -4.45, label = sprintf("Sensitivity cohort\n24/24 observed hours\nn = %s", format(fc$analyzed_ge24, big.mark = ",")), fill = "#ead3bc")
)

p <- ggplot() +
  coord_cartesian(xlim = c(-8.5, 8.6), ylim = c(-6.9, 14.4), clip = "off") +
  theme_void(base_family = "sans") +
  theme(
    plot.background = element_rect(fill = "#fbfaf6", color = NA),
    panel.background = element_rect(fill = "#fbfaf6", color = NA),
    plot.margin = margin(10, 36, 22, 36)
  )

for (i in seq_along(main_nodes)) {
  nd <- main_nodes[[i]]
  p <- p + flow_box(-2.75, 2.75, nd$y - 0.72, nd$y + 0.72, nd$label, nd$fill, text_size = 3.80, lineheight = 0.98)
}

for (i in seq_along(exclude_nodes)) {
  nd <- exclude_nodes[[i]]
  p <- p + flow_box(3.95, 8.40, nd$y - 0.72, nd$y + 0.72, nd$label, "#ffffff", size = 0.32, text_size = 3.36, lineheight = 0.95)
}

for (i in seq_along(branch_nodes)) {
  nd <- branch_nodes[[i]]
  p <- p + flow_box(nd$x - 2.35, nd$x + 2.35, nd$y - 0.70, nd$y + 0.70, nd$label, nd$fill, size = 0.32, text_size = 3.56, lineheight = 0.98)
}

for (y in c(13.6, 11.75, 9.9, 8.05, 6.2, 4.35, 2.5, 0.75)) {
  p <- p + annotate("segment", x = 0, xend = 0, y = y - 0.72, yend = y - 1.13,
                    linewidth = 0.42, color = "#3c4a52",
                    arrow = arrow(length = unit(0.14, "inches"), type = "closed"))
}

for (y in c(11.75, 9.9, 8.05, 6.2, 4.35, 2.5, 0.75, -1.50)) {
  p <- p + annotate("segment", x = 2.75, xend = 3.95, y = y, yend = y,
                    linewidth = 0.34, color = "#58656d",
                    arrow = arrow(length = unit(0.12, "inches"), type = "closed"))
}

p <- p +
  annotate("segment", x = 0, xend = -3.45, y = -2.22, yend = -3.68,
           linewidth = 0.34, color = "#55636b",
           arrow = arrow(length = unit(0.12, "inches"), type = "closed")) +
  annotate("segment", x = 0, xend = 3.45, y = -2.22, yend = -3.68,
           linewidth = 0.34, color = "#55636b",
           arrow = arrow(length = unit(0.12, "inches"), type = "closed")) +
  annotate("text", x = 0, y = -6.18,
           label = sprintf("Primary cohort required %s18/24 observed hourly invasive ABP-MAP values in the first 24 h; sensitivity cohorts used %s20/24 h and complete 24/24 h observation thresholds.", ge, ge),
           family = "sans", size = 3.32, color = "#33424a")

fig_png <- file.path(output_dir, "Figure_1_flow_diagram_major_revision.png")
fig_tiff <- file.path(output_dir, "Figure_1_flow_diagram_major_revision.tiff")
fig_pdf <- file.path(output_dir, "Figure_1_flow_diagram_major_revision.pdf")

ggsave(fig_png, p, width = 13.0, height = 12.4, dpi = 600, bg = "#fbfaf6")
ggsave(fig_tiff, p, width = 13.0, height = 12.4, dpi = 600, compression = "lzw", bg = "#fbfaf6")
ggsave(fig_pdf, p, width = 13.0, height = 12.4, bg = "#fbfaf6")

if (interactive()) {
  print(p)
}
}

fwrite(q4_profile_table, file.path(output_dir, "Table_R11_q4_deep_phenotype.csv"))
fwrite(q4_outcome_table, file.path(output_dir, "Table_R12_q4_survivors_vs_nonsurvivors.csv"))
fwrite(pheno_model_table, file.path(output_dir, "Table_R13_phenotype_augmented_model.csv"))
fwrite(pheno_model_overall, file.path(output_dir, "Table_R13_phenotype_augmented_model_overall.csv"))
fwrite(flow_counts, file.path(output_dir, "Table_R14_flow_counts_for_figure1.csv"))

try(dbDisconnect(con), silent = TRUE)

cat("\nQ4 deep phenotype analysis completed.\n")
if (build_figure1) cat("Figure 1 generation completed.\n")
cat("Output directory:", normalizePath(output_dir), "\n")
