# ============================================
# 01. Main analysis: cache generation and primary results
# Purpose:
# - Read the analysis cohort, hourly MAP data, and norepinephrine-equivalent exposure
# - Perform the 24-hour landmark analysis, multiple imputation, and primary Cox models
# - Generate the main cache and core outputs used for the manuscript
# ============================================

locate_current_script <- function() {
  frames <- sys.frames()
  if (length(frames) > 0) {
    for (i in rev(seq_along(frames))) {
      ofile <- frames[[i]]$ofile
      if (!is.null(ofile) && nzchar(ofile)) {
        return(normalizePath(ofile, winslash = "/", mustWork = FALSE))
      }
    }
  }

  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getSourceEditorContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) {
      return(normalizePath(ctx$path, winslash = "/", mustWork = FALSE))
    }
  }

  NA_character_
}

find_archive_root <- function(start_dir) {
  required_dirs <- c("01_SQL数据提取", "02_R统计复现", "03_图表", "04_文章", "05_附件")
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  repeat {
    if (all(dir.exists(file.path(current_dir, required_dirs)))) return(current_dir)
    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      stop("Unable to locate the archive root automatically. Run this script from the archive root or one of its subdirectories.")
    }
    current_dir <- parent_dir
  }
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
  library(tidyr)
  library(survival)
  library(mice)
  library(splines)
})
hour_cols   <- sprintf("h%02d", 0:23)
t0          <- 24
tau_icu     <- 30 * 24
m_imp       <- 5
maxit_imp   <- 10
seed_imp    <- 42

thr_hypo    <- 65
dt_hour     <- 1
output_dir <- "outputs_abp_map"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

show_png_in_plots <- function(path, title_text = NULL) {
  if (!interactive()) return(invisible(FALSE))
  if (!file.exists(path)) return(invisible(FALSE))
  if (!requireNamespace("png", quietly = TRUE) || !requireNamespace("grid", quietly = TRUE)) {
    return(invisible(FALSE))
  }

  img <- png::readPNG(path)
  grid::grid.newpage()
  grid::grid.raster(img)
  if (!is.null(title_text) && nzchar(title_text)) {
    grid::grid.text(
      label = title_text,
      x = 0.02, y = 0.98,
      just = c("left", "top"),
      gp = grid::gpar(fontsize = 12, fontface = "bold")
    )
  }
  invisible(TRUE)
}
safe_slug <- function(x) {
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^A-Za-z0-9_\\-]+", "", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (nchar(x) == 0) x <- "analysis"
  x
}

make_unique_prefix <- function(prefix, out_dir = ".") {
  pdf0 <- file.path(out_dir, paste0(prefix, ".pdf"))
  png0 <- file.path(out_dir, paste0(prefix, ".png"))
  if (!file.exists(pdf0) && !file.exists(png0)) return(file.path(out_dir, prefix))

  k <- 2
  repeat {
    pref_k <- paste0(prefix, "_v", k)
    pdf_k <- file.path(out_dir, paste0(pref_k, ".pdf"))
    png_k <- file.path(out_dir, paste0(pref_k, ".png"))
    if (!file.exists(pdf_k) && !file.exists(png_k)) return(file.path(out_dir, pref_k))
    k <- k + 1
  }
}
db <- "mimic4_v31"
con <- DBI::dbConnect(
  odbc::odbc(),
  dsn      = db,
  database = db,
  uid      = "postgres",
  pwd      = "postgres",
  server   = "localhost",
  port     = 5432
)
pop_sql <- "
SELECT *
FROM bdmcc.bdmcc_population;
"
pop_df <- as.data.table(DBI::dbGetQuery(con, pop_sql))
pop_df <- pop_df[order(stay_id)]
pop_df <- pop_df[, .SD[1], by = stay_id]

map_sql <- "
SELECT stay_id, hour, map_mean
FROM bdmcc.map_hourly_final_abp
WHERE hour BETWEEN 0 AND 23
ORDER BY stay_id, hour;
"
map_long <- as.data.table(DBI::dbGetQuery(con, map_sql))

ne_sql <- "
SELECT stay_id, ne_equiv_mean_0_24h
FROM bdmcc.vaso_ne_equiv_0_24h;
"
ne_df <- as.data.table(DBI::dbGetQuery(con, ne_sql))
map_wide <- dcast(map_long, stay_id ~ hour, value.var = "map_mean")
for (h in 0:23) {
  ch <- as.character(h)
  if (!ch %in% names(map_wide)) map_wide[, (ch) := NA_real_]
}
setcolorder(map_wide, c("stay_id", as.character(0:23)))
setnames(map_wide, old = as.character(0:23), new = hour_cols)
dat_abp <- pop_df %>%
  inner_join(map_wide, by = "stay_id") %>%
  left_join(ne_df, by = "stay_id") %>%
  as.data.table()

cat("ABP trajectory cohort size (before landmark):", nrow(dat_abp), "\n")
cat("Missing values in the NE variable:", sum(is.na(dat_abp$ne_equiv_mean_0_24h)), "\n")
stopifnot(all(c("day30_los","day30_outcome") %in% names(dat_abp)))

dat_lm <- dat_abp[day30_los > t0]
dat_lm[, time_lm_days := (pmin(day30_los, tau_icu) - t0) / 24]
dat_lm[, time_lm_days := pmax(time_lm_days, 1e-6)]
dat_lm[, event_lm := as.integer(day30_outcome == 1 & day30_los <= tau_icu & day30_los > t0)]

cat("ABP trajectory cohort size (after landmark):", nrow(dat_lm), "\n")
dat_lm[, n_hour_obs := rowSums(!is.na(as.matrix(.SD))), .SDcols = hour_cols]

cat("\nCohort size by hourly coverage threshold (after landmark):\n")
print(dat_lm[, .(
  n_ge_16 = sum(n_hour_obs >= 16),
  n_ge_18 = sum(n_hour_obs >= 18),
  n_ge_20 = sum(n_hour_obs >= 20),
  n_ge_24 = sum(n_hour_obs >= 24)
)])
has <- function(x, df) x[x %in% names(df)]

cov_m2_sofa <- unique(has(c(
  "age","gender","weight",
  "charlson","sofa",
  "itvtn_24h_vent_tag","itvtn_24h_rrt_tag",
  "drug_24h_sedative_tag"
), dat_lm))

cov_m2_sofa_ne <- unique(has(c(
  "age","gender","weight",
  "charlson","sofa",
  "itvtn_24h_vent_tag","itvtn_24h_rrt_tag",
  "drug_24h_sedative_tag",
  "ne_equiv_mean_0_24h"
), dat_lm))

cov_m2_saps <- unique(has(c(
  "age","gender","weight",
  "charlson","sapsii",
  "itvtn_24h_vent_tag","itvtn_24h_rrt_tag",
  "drug_24h_sedative_tag"
), dat_lm))

cat("\nModel 2 (SOFA) covariates:\n");  print(cov_m2_sofa)
cat("\nModel 2+ (SOFA + NE) covariates:\n"); print(cov_m2_sofa_ne)
cat("\nModel 2 (SAPSII) covariates:\n"); print(cov_m2_saps)
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
      if (length(idx) == 0) { dropped <- c(dropped, v); next }
      ref <- x[idx[1]]
      if (isTRUE(all(is.na(x) | x == ref))) { dropped <- c(dropped, v); next }
    }
  }

  cov_final <- setdiff(cov_use, dropped)
  list(covars = cov_final, dropped = unique(dropped))
}
calc_ttr_auc <- function(map_mat, thr = 65, dt = 1) {
  map_mat <- as.matrix(map_mat)
  if (is.null(dim(map_mat)) || ncol(map_mat) < 1) stop("map_mat must be a matrix with >=1 column.")
  if (anyNA(map_mat)) stop("map_mat contains NA; please impute or drop before computing TTR/AUC.")
  ttr <- rowMeans(map_mat < thr)
  gap <- thr - map_mat
  gap[gap < 0] <- 0
  auc <- rowSums(gap) * dt

  list(ttr = ttr, auc = auc)
}
run_mi_hours_only <- function(df, covars, hour_cols, m = 5, maxit = 10, seed = 42,
                              max_retry_drop_na = 2) {

  keep_cols <- unique(c("stay_id","time_lm_days","event_lm", covars, hour_cols))
  keep_cols <- keep_cols[keep_cols %in% names(df)]
  dat0 <- as.data.frame(df[, ..keep_cols])
  cat("\n[MI] Missing counts by hourly MAP column before imputation:\n")
  print(colSums(is.na(dat0[, hour_cols, drop = FALSE])))

  ini  <- mice(dat0, maxit = 0, printFlag = FALSE)
  meth <- ini$method
  pred <- ini$predictorMatrix
  meth[] <- ""
  meth[hour_cols] <- "pmm"
  for (v in c("stay_id","time_lm_days","event_lm")) if (v %in% names(meth)) meth[v] <- ""
  pred[,] <- 0
  cov_for_pred <- intersect(covars, colnames(pred))

  for (h in 0:23) {
    target <- sprintf("h%02d", h)
    if (!target %in% rownames(pred)) next

    nb <- (h-2):(h+2)
    nb <- nb[nb >= 0 & nb <= 23]
    nb_cols <- sprintf("h%02d", nb)
    nb_cols <- nb_cols[nb_cols %in% colnames(pred)]

    pred[target, nb_cols] <- 1
    pred[target, target]  <- 0

    if (length(cov_for_pred) > 0) pred[target, cov_for_pred] <- 1
    if ("time_lm_days" %in% colnames(pred)) pred[target, "time_lm_days"] <- 1
    if ("event_lm"     %in% colnames(pred)) pred[target, "event_lm"]     <- 1
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
    remove.constant  = TRUE,
    printFlag = TRUE
  )
  bad_ids <- c()
  for (i in 1:m) {
    di <- mice::complete(imp, i)
    na_row <- which(rowSums(is.na(di[, hour_cols, drop = FALSE])) > 0)
    if (length(na_row) > 0) {
      bad_ids <- union(bad_ids, di$stay_id[na_row])
    }
  }

  bad_ids <- unique(bad_ids)
  if (length(bad_ids) == 0) {
    return(list(imp = imp, dropped_stay_id = character(0)))
  }

  cat("\n[MI] WARNING: After imputation, some rows still have NA in hour_cols.\n")
  cat("[MI] Problematic stay_id count =", length(bad_ids), "\n")
  if (max_retry_drop_na <= 0) {
    cat("[MI] Reached max_retry_drop_na; returning current imp (may lead to NA downstream).\n")
    return(list(imp = imp, dropped_stay_id = bad_ids))
  }
  df2 <- as.data.table(df)
  df2 <- df2[!stay_id %in% bad_ids]

  cat("[MI] Dropping these stay_id and re-running MI. New N =", nrow(df2), "\n")

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
  return(res2)
}
extract_d1_fields <- function(d1_test) {
  cand <- c("statistic","df1","df2","dfcom","p.value","riv")
  out <- list()
  for (nm in cand) if (!is.null(d1_test[[nm]])) out[[nm]] <- d1_test[[nm]]

  if (length(out) >= 5 && !is.null(out[["p.value"]])) return(out)

  d1_out <- capture.output(print(d1_test))
  idx <- grep("~~", d1_out, fixed = TRUE)
  if (length(idx) == 0) idx <- grep("p", d1_out, ignore.case = TRUE)
  if (length(idx) == 0) stop("No parsable row was found in the D1 output.")

  d1_row <- d1_out[idx[1]]
  nums_chr <- regmatches(
    d1_row,
    gregexpr("[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?", d1_row, perl = TRUE)
  )[[1]]
  nums <- as.numeric(nums_chr)

  if (length(nums) >= 8) {
    list(statistic = nums[3], df1 = nums[4], df2 = nums[5], dfcom = nums[6], p.value = nums[7], riv = nums[8])
  } else {
    stop("Too few numeric values (<8) were parsed from the D1 output. Inspect d1_test:\n", d1_row)
  }
}
compute_exposure <- function(map_mat, type = c("mean_map","ttr_lt","auc_lt"), thr = 65, dt = 1) {
  type <- match.arg(type)
  if (type == "mean_map") {
    return(rowMeans(map_mat, na.rm = FALSE))
  }
  ba <- calc_ttr_auc(map_mat, thr = thr, dt = dt)
  if (type == "ttr_lt") return(ba$ttr)
  if (type == "auc_lt") return(ba$auc)
  stop("Unknown exposure type.")
}
fit_pooled_ns_and_D1_generic <- function(imp, covars, hour_cols,
                                         exposure_type = c("mean_map","ttr_lt","auc_lt"),
                                         thr = 65, dt = 1,
                                         time_var = "time_lm_days", event_var = "event_lm",
                                         df_spline = 4) {

  exposure_type <- match.arg(exposure_type)
  sc <- sanitize_covars_from_imp1(imp, covars)
  cov_use <- sc$covars
  vars_base <- unique(c(time_var, event_var, cov_use))

  x_term <- switch(exposure_type,
                   mean_map = "x_expo",
                   ttr_lt   = "x_expo",
                   auc_lt   = "x_expo")
  rhs_ns <- paste(c(sprintf("splines::ns(%s, df=%d)", x_term, df_spline), cov_use), collapse = " + ")
  fml_ns <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", rhs_ns))

  fit_ns <- with(imp, {
    map_mat <- do.call(cbind, mget(hour_cols, inherits = TRUE))
    x_expo  <- compute_exposure(map_mat, type = exposure_type, thr = thr, dt = dt)
    dd <- as.data.frame(mget(vars_base, inherits = TRUE))
    dd$x_expo <- x_expo
    coxph(fml_ns, data = dd, ties = "efron")
  })
  pooled_ns <- pool(fit_ns)
  rhs_lin <- paste(c(x_term, cov_use), collapse = " + ")
  fml_lin <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", rhs_lin))

  fit_lin <- with(imp, {
    map_mat <- do.call(cbind, mget(hour_cols, inherits = TRUE))
    x_expo  <- compute_exposure(map_mat, type = exposure_type, thr = thr, dt = dt)
    dd <- as.data.frame(mget(vars_base, inherits = TRUE))
    dd$x_expo <- x_expo
    coxph(fml_lin, data = dd, ties = "efron")
  })

  d1_test <- mice::D1(fit_ns, fit_lin)
  d1_parsed <- extract_d1_fields(d1_test)

  list(
    pooled_ns = pooled_ns,
    fit_ns = fit_ns,
    fit_lin = fit_lin,
    d1_test = d1_test,
    d1_parsed = d1_parsed,
    covars_used = cov_use,
    dropped = sc$dropped
  )
}
fit_pooled_quartile_map_mean <- function(imp, covars, hour_cols,
                                         time_var = "time_lm_days", event_var = "event_lm",
                                         cut_q_fixed = NULL) {

  sc <- sanitize_covars_from_imp1(imp, covars)
  cov_use <- sc$covars

  d1 <- complete(imp, 1)
  map_mat1 <- as.matrix(d1[, hour_cols, drop = FALSE])
  map24_mean1 <- rowMeans(map_mat1, na.rm = FALSE)
  if (is.null(cut_q_fixed)) {
  cut_q <- as.numeric(quantile(map24_mean1, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, type = 7))
} else {
  cut_q <- as.numeric(cut_q_fixed)
  if (length(cut_q) == 3) {
    cut_q <- c(min(map24_mean1, na.rm = TRUE), cut_q[1], cut_q[2], cut_q[3], max(map24_mean1, na.rm = TRUE))
  }
  if (length(cut_q) != 5) stop("cut_q_fixed must be length 5 (0/25/50/75/100%) or length 3 (25/50/75).")
}
  eps <- 1e-6
  for (i in 2:length(cut_q)) {
    if (!is.finite(cut_q[i])) cut_q[i] <- cut_q[i - 1] + eps
    if (cut_q[i] <= cut_q[i - 1]) cut_q[i] <- cut_q[i - 1] + eps
  }

  rhs_q <- paste(c("map_q", cov_use), collapse = " + ")
  fml_q <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", rhs_q))

  vars_base <- unique(c(time_var, event_var, cov_use))

  fit_q <- with(imp, {
    map_mat <- do.call(cbind, mget(hour_cols, inherits = TRUE))
    map24_mean <- rowMeans(map_mat, na.rm = FALSE)

    map_q <- cut(
      map24_mean,
      breaks = cut_q,
      include.lowest = TRUE,
      labels = c("Q1","Q2","Q3","Q4")
    )
    map_q <- relevel(map_q, ref = "Q2")

    dd <- as.data.frame(mget(vars_base, inherits = TRUE))
    dd$map_q <- map_q

    coxph(fml_q, data = dd, ties = "efron")
  })

  pooled_q <- pool(fit_q)
  list(cut_q = cut_q, pooled_q = pooled_q, covars_used = cov_use, dropped = sc$dropped)
}

rank_quartile_factor <- function(x) {
  r <- data.table::frank(x, ties.method = "first", na.last = "keep")

  br <- as.numeric(quantile(r, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, type = 1))
  eps <- 1e-6
  for (i in 2:length(br)) {
    if (!is.finite(br[i])) br[i] <- br[i - 1] + eps
    if (br[i] <= br[i - 1]) br[i] <- br[i - 1] + eps
  }

  q <- cut(r, breaks = br, include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4"))
  q <- relevel(q, ref = "Q2")
  q
}

fit_pooled_quartile_rank_generic <- function(imp, covars, hour_cols,
                                             exposure_type = c("ttr_lt","auc_lt"),
                                             thr = 65, dt = 1,
                                             time_var = "time_lm_days", event_var = "event_lm") {

  exposure_type <- match.arg(exposure_type)

  sc <- sanitize_covars_from_imp1(imp, covars)
  cov_use <- sc$covars

  rhs_q <- paste(c("q_expo", cov_use), collapse = " + ")
  fml_q <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", rhs_q))
  vars_base <- unique(c(time_var, event_var, cov_use))
  d1 <- complete(imp, 1)
  map_mat1 <- as.matrix(d1[, hour_cols, drop = FALSE])
  x1 <- compute_exposure(map_mat1, type = exposure_type, thr = thr, dt = dt)
  q1 <- rank_quartile_factor(x1)

  range_dt <- data.table(
    q = q1,
    x = x1
  )[, .(
    min_x = min(x, na.rm = TRUE),
    max_x = max(x, na.rm = TRUE)
  ), by = q][order(q)]

  fit_q <- with(imp, {
    map_mat <- do.call(cbind, mget(hour_cols, inherits = TRUE))
    x_expo  <- compute_exposure(map_mat, type = exposure_type, thr = thr, dt = dt)
    q_expo  <- rank_quartile_factor(x_expo)

    dd <- as.data.frame(mget(vars_base, inherits = TRUE))
    dd$q_expo <- q_expo
    coxph(fml_q, data = dd, ties = "efron")
  })

  pooled_q <- pool(fit_q)
  list(range_dt = range_dt, pooled_q = pooled_q, covars_used = cov_use, dropped = sc$dropped)
}
make_table_quartile_map <- function(pooled_q, cut_q) {
  sum_q <- summary(pooled_q, conf.int = TRUE, exponentiate = TRUE)
  sum_mapq <- sum_q[grepl("^map_q", sum_q$term), ]
  sum_mapq <- sum_mapq[match(c("map_qQ1","map_qQ3","map_qQ4"), sum_mapq$term), ]

  q1_range <- sprintf("[%.2f, %.2f)", cut_q[1], cut_q[2])
  q3_range <- sprintf("[%.2f, %.2f)", cut_q[3], cut_q[4])
  q4_range <- sprintf("[%.2f, %.2f]", cut_q[4], cut_q[5])

  table2 <- data.frame(
    Exposure = c("MAP Q1 vs Q2", "MAP Q3 vs Q2", "MAP Q4 vs Q2"),
    MAP_Quartile_Range_mmHg = c(q1_range, q3_range, q4_range),
    HR  = sum_mapq$estimate,
    LCI = sum_mapq$`2.5 %`,
    UCI = sum_mapq$`97.5 %`,
    P   = sum_mapq$p.value,
    stringsAsFactors = FALSE
  )
  table2$`HR (95% CI)` <- sprintf("%.2f (%.2f–%.2f)", table2$HR, table2$LCI, table2$UCI)
  table2[, c("Exposure","MAP_Quartile_Range_mmHg","HR (95% CI)","P")]
}

make_table_quartile_rank <- function(pooled_q, range_dt, label_prefix,
                                     unit_label = "") {
  sum_q <- summary(pooled_q, conf.int = TRUE, exponentiate = TRUE)
  sum_qx <- sum_q[grepl("^q_expo", sum_q$term), ]
  sum_qx <- sum_qx[match(c("q_expoQ1","q_expoQ3","q_expoQ4"), sum_qx$term), ]
  rng <- as.data.frame(range_dt)
  rng <- rng[match(c("Q1","Q3","Q4"), rng$q), ]

  rng_txt <- if (unit_label == "") {
    sprintf("[%.4f, %.4f]", rng$min_x, rng$max_x)
  } else {
    sprintf("[%.4f, %.4f] %s", rng$min_x, rng$max_x, unit_label)
  }

  table2 <- data.frame(
    Exposure = c(paste0(label_prefix, " Q1 vs Q2"),
                 paste0(label_prefix, " Q3 vs Q2"),
                 paste0(label_prefix, " Q4 vs Q2")),
    Quartile_Range = rng_txt,
    HR  = sum_qx$estimate,
    LCI = sum_qx$`2.5 %`,
    UCI = sum_qx$`97.5 %`,
    P   = sum_qx$p.value,
    stringsAsFactors = FALSE
  )
  table2$`HR (95% CI)` <- sprintf("%.2f (%.2f–%.2f)", table2$HR, table2$LCI, table2$UCI)
  table2[, c("Exposure","Quartile_Range","HR (95% CI)","P")]
}
make_table_linear <- function(pooled_lin, scale = 1, label = "Exposure") {
  sum_lin <- summary(pooled_lin, conf.int = TRUE, exponentiate = TRUE)
  row <- sum_lin[sum_lin$term == "x_expo", ]
  if (nrow(row) != 1) stop("Cannot find term x_expo in pooled linear model.")
  est <- log(row$estimate) * scale
  se  <- ( (log(row$`97.5 %`) - log(row$estimate)) / 1.96 ) * scale
  hr  <- exp(est)
  lci <- exp(est - 1.96*se)
  uci <- exp(est + 1.96*se)

  out <- data.frame(
    Exposure = label,
    `HR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", hr, lci, uci),
    P = row$p.value,
    stringsAsFactors = FALSE
  )
  out
}
plot_spline_hr_generic <- function(imp, covars, hour_cols,
                                   exposure_type = c("mean_map","ttr_lt","auc_lt"),
                                   thr = 65, dt = 1,
                                   time_var = "time_lm_days", event_var = "event_lm",
                                   df_spline = 4, grid_n = 200,
                                   q_trim = c(0.01, 0.99),
                                   xlab = NULL, main = NULL,
                                   out_prefix = "Figure_Spline_HR",
                                   out_dir = ".") {

  exposure_type <- match.arg(exposure_type)

  d_fig <- mice::complete(imp, 1)
  map_mat <- as.matrix(d_fig[, hour_cols, drop = FALSE])

  x_expo <- compute_exposure(map_mat, type = exposure_type, thr = thr, dt = dt)
  d_fig$x_expo <- x_expo
  for (v in covars) {
    if (!v %in% names(d_fig)) next
    if (is.character(d_fig[[v]])) d_fig[[v]] <- factor(d_fig[[v]])
  }
  drop_single_level <- c()
  drop_constant_num <- c()
  for (v in covars) {
    if (!v %in% names(d_fig)) next
    x <- d_fig[[v]]
    if (is.factor(x) && nlevels(x) < 2) drop_single_level <- c(drop_single_level, v)
    if (is.numeric(x) || is.integer(x)) {
      idx <- which(!is.na(x))
      if (length(idx) == 0) { drop_constant_num <- c(drop_constant_num, v) } else {
        ref <- x[idx[1]]
        if (isTRUE(all(is.na(x) | x == ref))) drop_constant_num <- c(drop_constant_num, v)
      }
    }
  }
  cov_use <- setdiff(covars, unique(c(drop_single_level, drop_constant_num)))

  if (length(drop_single_level) > 0) message("[Figure] Dropped single-level factor covariate(s): ", paste(unique(drop_single_level), collapse = ", "))
  if (length(drop_constant_num) > 0) message("[Figure] Dropped constant numeric covariate(s): ", paste(unique(drop_constant_num), collapse = ", "))
  rhs <- paste(c(sprintf("splines::ns(x_expo, df=%d)", df_spline), cov_use), collapse = " + ")
  fml <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", rhs))
  fit <- coxph(fml, data = d_fig, ties = "efron")
  make_typical_row <- function(df, covars) {
    out <- list()
    for (v in covars) {
      if (!v %in% names(df)) next
      x <- df[[v]]
      if (is.factor(x)) {
        tb <- table(x)
        mode_level <- names(tb)[which.max(tb)]
        out[[v]] <- factor(mode_level, levels = levels(x))
      } else if (is.numeric(x) || is.integer(x)) {
        out[[v]] <- stats::median(x, na.rm = TRUE)
      } else if (is.logical(x)) {
        out[[v]] <- FALSE
      } else {
        out[[v]] <- stats::median(as.numeric(x), na.rm = TRUE)
      }
    }
    as.data.frame(out, stringsAsFactors = FALSE)
  }
  typical <- make_typical_row(d_fig, cov_use)
  x_ref <- as.numeric(stats::median(d_fig$x_expo, na.rm = TRUE))
  x_lo  <- as.numeric(stats::quantile(d_fig$x_expo, q_trim[1], na.rm = TRUE))
  x_hi  <- as.numeric(stats::quantile(d_fig$x_expo, q_trim[2], na.rm = TRUE))
  x_grid <- seq(x_lo, x_hi, length.out = grid_n)

  nd_grid <- typical[rep(1, length(x_grid)), , drop = FALSE]
  nd_grid$x_expo <- x_grid
  nd_ref <- typical
  nd_ref$x_expo <- x_ref
  for (v in cov_use) {
    if (!v %in% names(d_fig)) next
    if (is.factor(d_fig[[v]])) {
      nd_grid[[v]] <- factor(nd_grid[[v]], levels = levels(d_fig[[v]]))
      nd_ref[[v]]  <- factor(nd_ref[[v]],  levels = levels(d_fig[[v]]))
    }
  }
  Xg <- model.matrix(delete.response(terms(fit)), data = nd_grid)
  Xr <- model.matrix(delete.response(terms(fit)), data = nd_ref)

  beta <- coef(fit)
  V <- vcov(fit)

  common <- intersect(colnames(Xg), names(beta))
  Xg <- Xg[, common, drop = FALSE]
  Xr <- Xr[, common, drop = FALSE]
  beta <- beta[common]
  V <- V[common, common, drop = FALSE]

  D <- Xg - matrix(Xr[1, ], nrow = nrow(Xg), ncol = ncol(Xg), byrow = TRUE)
  lp <- as.vector(D %*% beta)
  se <- sqrt(pmax(0, rowSums((D %*% V) * D)))

  hr  <- exp(lp)
  lcl <- exp(lp - 1.96 * se)
  ucl <- exp(lp + 1.96 * se)

  df_plot <- data.frame(X = x_grid, HR = hr, LCL = lcl, UCL = ucl)

  if (is.null(xlab)) xlab <- "Exposure"
  if (is.null(main)) main <- ""
  if (!grepl(.Platform$file.sep, out_prefix, fixed = TRUE)) {
    out_prefix <- file.path(out_dir, out_prefix)
  }
  out_dir2 <- dirname(out_prefix)
  base2 <- basename(out_prefix)
  out_prefix_unique <- make_unique_prefix(base2, out_dir = out_dir2)

  pdf_file <- paste0(out_prefix_unique, ".pdf")
  png_file <- paste0(out_prefix_unique, ".png")
  tiff_file <- paste0(out_prefix_unique, ".tiff")
  cairo_pdf(pdf_file, width = 8, height = 6, family = "Arial")
  par(mar = c(5, 5, 2, 2), mgp = c(3.5, 1, 0), las = 1,
      cex.lab = 1.3, cex.axis = 1.1, font.lab = 2)

  plot(df_plot$X, df_plot$HR, type = "l", lwd = 2.5, col = "#3C5488",
       xlab = xlab, ylab = "Hazard ratio (95% CI)",
       ylim = range(c(df_plot$LCL, df_plot$UCL), na.rm = TRUE),
       bty = "l", axes = FALSE)

  polygon(c(df_plot$X, rev(df_plot$X)),
          c(df_plot$LCL, rev(df_plot$UCL)),
          col = rgb(77, 187, 213, alpha = 50, maxColorValue = 255),
          border = NA)

  lines(df_plot$X, df_plot$HR, lwd = 2.5, col = "#3C5488")
  abline(h = 1, lty = 2, lwd = 1.5, col = "#E64B35")
  abline(v = x_ref, lty = 3, lwd = 1.5, col = "#00A087")

  axis(1, lwd = 1.5, col = "black")
  axis(2, lwd = 1.5, col = "black")
  box(lwd = 1.5)

  dev.off()
  png(png_file, width = 8, height = 6, units = "in", res = 600, type = "cairo")
  par(mar = c(5, 5, 2, 2), mgp = c(3.5, 1, 0), las = 1,
      cex.lab = 1.3, cex.axis = 1.1, font.lab = 2)

  plot(df_plot$X, df_plot$HR, type = "l", lwd = 2.5, col = "#3C5488",
       xlab = xlab, ylab = "Hazard ratio (95% CI)",
       ylim = range(c(df_plot$LCL, df_plot$UCL), na.rm = TRUE),
       bty = "l", axes = FALSE)

  polygon(c(df_plot$X, rev(df_plot$X)),
          c(df_plot$LCL, rev(df_plot$UCL)),
          col = rgb(77, 187, 213, alpha = 50, maxColorValue = 255),
          border = NA)

  lines(df_plot$X, df_plot$HR, lwd = 2.5, col = "#3C5488")
  abline(h = 1, lty = 2, lwd = 1.5, col = "#E64B35")
  abline(v = x_ref, lty = 3, lwd = 1.5, col = "#00A087")

  axis(1, lwd = 1.5, col = "black")
  axis(2, lwd = 1.5, col = "black")
  box(lwd = 1.5)

  dev.off()
  tiff(tiff_file, width = 8, height = 6, units = "in", res = 600, compression = "lzw")
  par(mar = c(5, 5, 2, 2), mgp = c(3.5, 1, 0), las = 1,
      cex.lab = 1.3, cex.axis = 1.1, font.lab = 2)

  plot(df_plot$X, df_plot$HR, type = "l", lwd = 2.5, col = "#3C5488",
       xlab = xlab, ylab = "Hazard ratio (95% CI)",
       ylim = range(c(df_plot$LCL, df_plot$UCL), na.rm = TRUE),
       bty = "l", axes = FALSE)

  polygon(c(df_plot$X, rev(df_plot$X)),
          c(df_plot$LCL, rev(df_plot$UCL)),
          col = rgb(77, 187, 213, alpha = 50, maxColorValue = 255),
          border = NA)

  lines(df_plot$X, df_plot$HR, lwd = 2.5, col = "#3C5488")
  abline(h = 1, lty = 2, lwd = 1.5, col = "#E64B35")
  abline(v = x_ref, lty = 3, lwd = 1.5, col = "#00A087")

  axis(1, lwd = 1.5, col = "black")
  axis(2, lwd = 1.5, col = "black")
  box(lwd = 1.5)

  dev.off()

  cat("\n[Figure] Publication-quality figures saved:\n")
  cat("  PDF (vector): ", pdf_file, "\n")
  cat("  PNG (600 DPI): ", png_file, "\n")
  cat("  TIFF (600 DPI): ", tiff_file, "\n", sep = "")

  invisible(list(
    df_plot = df_plot,
    fit = fit,
    x_ref = x_ref,
    covars_used = cov_use,
    pdf = pdf_file,
    png = png_file,
    tiff = tiff_file
  ))
}
run_analysis_bundle <- function(df_lm, hours_threshold, covars, hour_cols, label,
                                m = 5, maxit = 10, seed = 42,
                                thr_hypo = 65, dt = 1,
                                do_ttr = TRUE, do_ttr_quartile = FALSE, do_plots = TRUE,
                                out_dir = ".",
                                map_cut_q = NULL) {
  df_raw <- df_lm[n_hour_obs >= hours_threshold]
  cat("\n[Bundle] ", label, " | N_raw=", nrow(df_raw),
      " events=", sum(df_raw$event_lm, na.rm = TRUE), "\n", sep="")
  mi_res <- run_mi_hours_only(df_raw, covars, hour_cols, m = m, maxit = maxit, seed = seed,
                              max_retry_drop_na = 2)
  imp <- mi_res$imp
  drop_ids <- mi_res$dropped_stay_id
  
  df_used <- df_raw
  if (length(drop_ids) > 0) {
    df_used <- df_raw[!stay_id %in% drop_ids]
  }
  
  cat("[Bundle] N_used(after MI NA-check)=", nrow(df_used),
      " dropped_for_MI=", length(drop_ids), "\n", sep="")
  res_q_map <- fit_pooled_quartile_map_mean(imp, covars, hour_cols, cut_q_fixed = map_cut_q)
  table2_map <- make_table_quartile_map(res_q_map$pooled_q, res_q_map$cut_q)
  res_ns_map <- fit_pooled_ns_and_D1_generic(
    imp = imp, 
    covars = covars, 
    hour_cols = hour_cols,
    exposure_type = "mean_map",
    thr = thr_hypo, 
    dt = dt,
    time_var = "time_lm_days", 
    event_var = "event_lm",
    df_spline = 4
  )
  res_ns_ttr <- NULL
  table_ttr_lin <- data.frame()
  table2_ttr_q <- data.frame()
  
  if (isTRUE(do_ttr)) {
    res_ns_ttr <- fit_pooled_ns_and_D1_generic(
      imp = imp,
      covars = covars,
      hour_cols = hour_cols,
      exposure_type = "ttr_lt",
      thr = thr_hypo, 
      dt = dt,
      time_var = "time_lm_days", 
      event_var = "event_lm",
      df_spline = 4
    )
    pooled_ttr_lin <- pool(res_ns_ttr$fit_lin)
    table_ttr_lin <- make_table_linear(
      pooled_ttr_lin, 
      scale = 0.10, 
      label = "TTR<65 (per 10%)"
    )
    if (isTRUE(do_ttr_quartile)) {
      res_q_ttr_list <- fit_pooled_quartile_rank_generic(
        imp = imp, 
        covars = covars, 
        hour_cols = hour_cols,
        exposure_type = "ttr_lt",
        thr = thr_hypo, 
        dt = dt,
        time_var = "time_lm_days", 
        event_var = "event_lm"
      )
      table2_ttr_q <- make_table_quartile_rank(
        res_q_ttr_list$pooled_q, 
        res_q_ttr_list$range_dt, 
        label_prefix = "TTR<65"
      )
    }
  }
  fig_map <- fig_ttr <- NULL
  if (isTRUE(do_plots)) {
    fig_map_prefix <- file.path(out_dir, paste0("Figure_Spline_HR_vs_meanMAP_", safe_slug(label)))
    fig_map <- plot_spline_hr_generic(
      imp = imp,
      covars = res_ns_map$covars_used,
      hour_cols = hour_cols,
      exposure_type = "mean_map",
      thr = thr_hypo, dt = dt,
      xlab = "TWA-MAP (mmHg) during first 24 hours",
      main = "Natural spline (df=4): HR vs TWA-MAP (ref=median)",
      out_prefix = fig_map_prefix,
      out_dir = out_dir
    )
    
    if (isTRUE(do_ttr)) {
      fig_ttr_prefix <- file.path(out_dir, paste0("Figure_Spline_HR_vs_TTRlt65_", safe_slug(label)))
      fig_ttr <- plot_spline_hr_generic(
        imp = imp,
        covars = res_ns_ttr$covars_used,
        hour_cols = hour_cols,
        exposure_type = "ttr_lt",
        thr = thr_hypo, dt = dt,
        xlab = "TTR<65 (proportion of hours with MAP<65)",
        main = "Natural spline (df=4): HR vs TTR<65 (ref=median)",
        out_prefix = fig_ttr_prefix,
        out_dir = out_dir
      )
    }
  }
  write.csv(table2_map, file.path(out_dir, paste0("Table2_MAP_", safe_slug(label), ".csv")), row.names = FALSE)
  if (isTRUE(do_ttr)) {
    write.csv(table_ttr_lin, file.path(out_dir, paste0("Table2_TTR_linear_", safe_slug(label), ".csv")), row.names = FALSE)
    if (isTRUE(do_ttr_quartile)) {
      write.csv(table2_ttr_q, file.path(out_dir, paste0("Table2_TTR_quartile_", safe_slug(label), ".csv")), row.names = FALSE)
    }
  }
  get_cell_chr <- function(df, i, col) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) < i || !(col %in% names(df))) return(NA_character_)
    as.character(df[[col]][i])
  }
  get_cell_num <- function(df, i, col) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) < i || !(col %in% names(df))) return(NA_real_)
    suppressWarnings(as.numeric(df[[col]][i]))
  }
  
  sum_row <- data.frame(
    Analysis = label,
    N_raw = nrow(df_raw),
    N_used = nrow(df_used),
    Dropped_for_MI = length(drop_ids),
    Events = sum(df_used$event_lm, na.rm = TRUE),
    MAP_Q1_vs_Q2 = get_cell_chr(table2_map, 1, "HR (95% CI)"),
    MAP_P_Q1     = get_cell_num(table2_map, 1, "P"),
    MAP_Q3_vs_Q2 = get_cell_chr(table2_map, 2, "HR (95% CI)"),
    MAP_P_Q3     = get_cell_num(table2_map, 2, "P"),
    MAP_Q4_vs_Q2 = get_cell_chr(table2_map, 3, "HR (95% CI)"),
    MAP_P_Q4     = get_cell_num(table2_map, 3, "P"),
    MAP_Nonlinearity_P = if (!is.null(res_ns_map$d1_parsed$p.value)) res_ns_map$d1_parsed$p.value else NA_real_,
    TTR_per10pct_HR = if (isTRUE(do_ttr)) get_cell_chr(table_ttr_lin, 1, "HR (95% CI)") else NA_character_,
    TTR_per10pct_P  = if (isTRUE(do_ttr)) get_cell_num(table_ttr_lin, 1, "P") else NA_real_,
    TTR_Nonlinearity_P = if (isTRUE(do_ttr) && !is.null(res_ns_ttr$d1_parsed$p.value)) res_ns_ttr$d1_parsed$p.value else NA_real_,
    
    stringsAsFactors = FALSE
  )
  
  return(list(
    imp = imp,
    dropped_ids = drop_ids,
    full_data = mice::complete(imp, 1), 
    res_q_map = res_q_map,      
    table2_map = table2_map,
    res_ns_map = res_ns_map,
    table_ttr_lin = table_ttr_lin,
    table2_ttr_q = table2_ttr_q,
    res_ns_ttr = res_ns_ttr,
    fig_map = fig_map,
    fig_ttr = fig_ttr,
    summary_row = sum_row
  ))
}

cat("\n[Step 19] Starting Main and Sensitivity Analyses...\n")
main_res <- run_analysis_bundle(
  df_lm = dat_lm,
  hours_threshold = 18,
  covars = cov_m2_sofa,
  hour_cols = hour_cols,
  label = "Main_ge18h_Model2_SOFA",
  m = m_imp, maxit = maxit_imp, seed = seed_imp,
  thr_hypo = thr_hypo, dt = dt_hour,
  do_ttr = TRUE,
  do_ttr_quartile = FALSE,
  do_plots = TRUE,
  out_dir = output_dir
)
main_cut_q <- main_res$res_q_map$cut_q
cat("[Info] Fixed MAP quartile cut-points (main cohort) = ",
    paste(sprintf("%.2f", main_cut_q), collapse = ", "), "\n", sep="")
main_res_ne <- run_analysis_bundle(
  df_lm = dat_lm,
  hours_threshold = 18,
  covars = cov_m2_sofa_ne,
  hour_cols = hour_cols,
  label = "Main_ge18h_Model2plus_SOFA_NE",
  m = m_imp, maxit = maxit_imp, seed = seed_imp,
  thr_hypo = thr_hypo, dt = dt_hour,
  do_ttr = TRUE,
  do_ttr_quartile = FALSE,
  do_plots = TRUE,
  out_dir = output_dir,
  map_cut_q = main_cut_q
)
sens20_res <- run_analysis_bundle(
  df_lm = dat_lm,
  hours_threshold = 20,
  covars = cov_m2_sofa,
  hour_cols = hour_cols,
  label = "Sensitivity_ge20h_Model2_SOFA",
  m = m_imp, maxit = maxit_imp, seed = seed_imp,
  thr_hypo = thr_hypo, dt = dt_hour,
  do_ttr = FALSE,
  do_plots = TRUE,
  out_dir = output_dir,
  map_cut_q = main_cut_q
)
sens24_res <- run_analysis_bundle(
  df_lm = dat_lm,
  hours_threshold = 24,
  covars = cov_m2_sofa,
  hour_cols = hour_cols,
  label = "Sensitivity_ge24h_Model2_SOFA",
  m = m_imp, maxit = maxit_imp, seed = seed_imp,
  thr_hypo = thr_hypo, dt = dt_hour,
  do_ttr = FALSE,
  do_plots = TRUE,
  out_dir = output_dir,
  map_cut_q = main_cut_q
)
sensSAPS_res <- run_analysis_bundle(
  df_lm = dat_lm,
  hours_threshold = 18,
  covars = cov_m2_saps,
  hour_cols = hour_cols,
  label = "Sensitivity_SAPSII_replacing_SOFA_ge18h",
  m = m_imp, maxit = maxit_imp, seed = seed_imp,
  thr_hypo = thr_hypo, dt = dt_hour,
  do_ttr = FALSE,
  do_plots = TRUE,
  out_dir = output_dir,
  map_cut_q = main_cut_q
)
final_data_csv <- file.path(output_dir, "final_imputed_data_for_table1.csv")
write.csv(main_res$full_data, final_data_csv, row.names = FALSE)
cat("\n[Success] Imputed Table 1 source data exported to: ", final_data_csv, "\n")

cat("\n[Step 20] Summarizing all results into Table S...\n")

tableS_all <- data.table::rbindlist(
  list(
    main_res$summary_row,
    main_res_ne$summary_row,
    sens20_res$summary_row,
    sens24_res$summary_row,
    sensSAPS_res$summary_row
  ),
  fill = TRUE
)

tableS_csv <- file.path(output_dir, "TableS_summary_MAP_TTR.csv")
write.csv(tableS_all, tableS_csv, row.names = FALSE)

cat("[Output] TableS summary saved: ", normalizePath(tableS_csv), "\n")
cat("\n[Step 21] Exporting cache for table1.R...\n")

cache_dir <- file.path(output_dir, "cache")
if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

cache_data <- list(
  dat_lm18 = dat_lm[n_hour_obs >= 18],
  imp1_hours = mice::complete(main_res$imp, 1), 
  hour_cols = hour_cols,
  cut_q_main = main_res$res_q_map$cut_q  
)

saveRDS(cache_data, file.path(cache_dir, "cache_main_ge18h_Model2_SOFA_minimal.rds"))
cat("[Success] Cache file generated with cut_q_main length:", length(cache_data$cut_q_main), "\n")

cat("\n===== FINAL N & EVENTS CHECK =====\n")
print(tableS_all[, .(Analysis, N_raw, N_used, Dropped_for_MI, Events)])

if (interactive()) {
  preview_pngs <- sort(list.files(output_dir, pattern = "^Figure_.*\\.png$", full.names = TRUE))
  preview_pngs <- preview_pngs[file.exists(preview_pngs)]
  if (length(preview_pngs) > 0) {
    cat("\n[Preview] Displaying generated figures in the RStudio Plots pane...\n")
    invisible(lapply(preview_pngs, function(fig) {
      show_png_in_plots(fig, basename(fig))
    }))
  }
}

try(DBI::dbDisconnect(con), silent = TRUE)

cat("\n[Done] All analyses completed. You can now run table1.R safely.\n")
