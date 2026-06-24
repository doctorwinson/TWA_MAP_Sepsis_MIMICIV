# ============================================
# 02. Table 1 and baseline tables
# Purpose:
# - Prepare the baseline dataset from the main analysis cache
# - Generate the manuscript Table 1 dataset and the full baseline table with SMDs
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
  library(data.table)
  library(dplyr)
  library(purrr)
  library(strong)
  library(DBI)
  library(odbc)
})
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

hour_cols <- sprintf("h%02d", 0:23)
cache_file <- "outputs_abp_map/cache/cache_main_ge18h_Model2_SOFA_minimal.rds"
stopifnot(file.exists(cache_file))

cache <- readRDS(cache_file)
stopifnot(!is.null(cache$dat_lm18))
dat_lm18 <- as.data.table(cache$dat_lm18)
stopifnot(!is.null(cache$imp1_hours))
imp1_hours <- as.data.table(cache$imp1_hours)
hour_cols_cache <- cache$hour_cols
cut_q_main <- as.numeric(cache$cut_q_main)

stopifnot(length(hour_cols_cache) == 24)
stopifnot(all(hour_cols_cache %in% names(imp1_hours)))
stopifnot(all(c("stay_id", hour_cols_cache) %in% names(imp1_hours)))
stopifnot(all("stay_id" %in% names(dat_lm18)))
stopifnot(length(cut_q_main) == 5)
imp1_hours[, map24_mean := rowMeans(.SD, na.rm = FALSE), .SDcols = hour_cols_cache]
na_mapmean <- sum(is.na(imp1_hours$map24_mean))
cat("[Check] imp#1 map24_mean NA =", na_mapmean, "\n")

imp1_hours[, MAP_group := cut(
  map24_mean,
  breaks = cut_q_main,
  include.lowest = TRUE,
  labels = c("Q1", "Q2", "Q3", "Q4")
)]

imp1_hours[, MAP_group := factor(MAP_group, levels = c("Q1","Q2","Q3","Q4"))]
dat_lm18 <- merge(
  dat_lm18,
  imp1_hours[, .(stay_id, map24_mean, MAP_group)],
  by = "stay_id",
  all.x = TRUE
)
dat_lm18 <- dat_lm18[!is.na(MAP_group)]
na_grp <- sum(is.na(dat_lm18$MAP_group))
cat("[Check] dat_lm18 N =", nrow(dat_lm18),
    " Events =", if ("event_lm" %in% names(dat_lm18)) sum(dat_lm18$event_lm, na.rm = TRUE) else NA,
    " NA_group =", na_grp, "\n")
print(table(dat_lm18$MAP_group, useNA = "ifany"))

stopifnot(na_grp == 0)
vaso_query <- "SELECT * FROM bdmcc.vaso_ne_equiv_0_24h"
vaso_data <- setDT(dbGetQuery(con, vaso_query))
ne_vars_to_remove <- c("ne_equiv_auc_0_24h", "ne_equiv_mean_0_24h", "ne_equiv_max_0_24h")
existing_ne_vars <- intersect(ne_vars_to_remove, names(dat_lm18))
if (length(existing_ne_vars) > 0) {
  cat("[Warning] Removing existing NE variables from dat_lm18:", paste(existing_ne_vars, collapse=", "), "\n")
  dat_lm18[, (existing_ne_vars) := NULL]
}

dat_lm18 <- merge(
  dat_lm18,
  vaso_data[, .(stay_id, ne_equiv_auc_0_24h, ne_equiv_mean_0_24h, ne_equiv_max_0_24h)],
  by = "stay_id",
  all.x = TRUE
)
keep_vars <- intersect(c(
  "stay_id", "MAP_group",
  "age", "gender", "weight",
  "vs_24h_temp_first", "vs_24h_heart_rate_first", "vs_24h_resp_rate_first", "map24_mean",
  "lab_fst_24h_wbc_first", "lab_fst_24h_hemoglobin_first", "lab_fst_24h_hct_first", "lab_fst_24h_platelet_first",
  "lab_fst_24h_bun_first", "lab_fst_24h_creatinine_first", "lab_fst_24h_k_first", "lab_fst_24h_na_first",
  "lab_fst_24h_cl_first", "lab_fst_24h_pt_first","lab_fst_24h_lactate_first",
  "charlson", "sofa", "sapsii",
  "icd_hypertension", "icd_diabetes", "icd_hf", "icd_copd", "icd_stroke",
  "icd_renal", "icd_lc", "icd_cad",
  "itvtn_24h_vent_tag", "itvtn_24h_rrt_tag", "drug_24h_sedative_tag",
  "ne_equiv_auc_0_24h", "ne_equiv_mean_0_24h", "ne_equiv_max_0_24h",
  "time_lm_days", "event_lm"
), names(dat_lm18))

table1_df <- dat_lm18[, ..keep_vars]

cat("\n[Table1] Keep vars:", length(keep_vars), "\n")
cat("[Table1] N =", nrow(table1_df),
    " Events =", if ("event_lm" %in% names(table1_df)) sum(table1_df$event_lm, na.rm = TRUE) else NA, "\n")
print(table(table1_df$MAP_group, useNA = "ifany"))
col_name_gp <- "Group"
gp_name <- c("Q1", "Q2", "Q3", "Q4")

header_name <- c(
  "Demographic Characteristics",
  "Vital Signs",
  "Laboratory Tests",
  "Disease Severity Scoring Systems",
  "Comorbid Conditions",
  "Interventions",
  "Vasopressor Exposure (0-24h)"
)
header_loc <- c(1, 5, 10, 22, 26, 35, 38)

df0 <- as.data.frame(table1_df)
cat("\n[Diagnostic] NE variables in df0:\n")
cat("  Column includes ne_equiv_auc_0_24h:", "ne_equiv_auc_0_24h" %in% names(df0), "\n")
cat("  Column includes ne_equiv_mean_0_24h:", "ne_equiv_mean_0_24h" %in% names(df0), "\n")
cat("  Column includes ne_equiv_max_0_24h:", "ne_equiv_max_0_24h" %in% names(df0), "\n")
if ("ne_equiv_mean_0_24h" %in% names(df0)) {
  cat("  Missing values in ne_equiv_mean_0_24h:", sum(is.na(df0$ne_equiv_mean_0_24h)), "/", nrow(df0), "\n")
  cat("  Range of ne_equiv_mean_0_24h:", paste(range(df0$ne_equiv_mean_0_24h, na.rm=TRUE), collapse=" ~ "), "\n")
}
safe_factor01 <- function(x) {
  factor(x, levels = c(1, 0), labels = c("YES", "NO"))
}

df <- df0 %>%
  mutate(
    MAP_group = factor(MAP_group, levels = c("Q1","Q2","Q3","Q4"), labels = gp_name)
  ) %>%
  { if ("gender" %in% names(.)) mutate(., gender = factor(gender, levels = c("F","M"), labels = c("Female","Male"))) else . } %>%
  { if ("itvtn_24h_vent_tag" %in% names(.)) mutate(., itvtn_24h_vent_tag = safe_factor01(itvtn_24h_vent_tag)) else . } %>%
  { if ("itvtn_24h_rrt_tag" %in% names(.)) mutate(., itvtn_24h_rrt_tag  = safe_factor01(itvtn_24h_rrt_tag)) else . } %>%
  { if ("drug_24h_sedative_tag" %in% names(.)) mutate(., drug_24h_sedative_tag = safe_factor01(drug_24h_sedative_tag)) else . } %>%
  { if ("icd_hypertension" %in% names(.)) mutate(., icd_hypertension = safe_factor01(icd_hypertension)) else . } %>%
  { if ("icd_diabetes" %in% names(.)) mutate(., icd_diabetes     = safe_factor01(icd_diabetes)) else . } %>%
  { if ("icd_hf" %in% names(.)) mutate(., icd_hf           = safe_factor01(icd_hf)) else . } %>%
  { if ("icd_copd" %in% names(.)) mutate(., icd_copd         = safe_factor01(icd_copd)) else . } %>%
  { if ("icd_stroke" %in% names(.)) mutate(., icd_stroke       = safe_factor01(icd_stroke)) else . } %>%
  { if ("icd_renal" %in% names(.)) mutate(., icd_renal        = safe_factor01(icd_renal)) else . } %>%
  { if ("icd_lc" %in% names(.)) mutate(., icd_lc           = safe_factor01(icd_lc)) else . } %>%
  { if ("icd_cad" %in% names(.)) mutate(., icd_cad          = safe_factor01(icd_cad)) else . } %>%
  { if ("ne_equiv_auc_0_24h" %in% names(.)) . else mutate(., ne_equiv_auc_0_24h = NA_real_) } %>%
  { if ("ne_equiv_mean_0_24h" %in% names(.)) . else mutate(., ne_equiv_mean_0_24h = NA_real_) } %>%
  { if ("ne_equiv_max_0_24h" %in% names(.)) . else mutate(., ne_equiv_max_0_24h = NA_real_) } %>%
  dplyr::rename(
    Group = MAP_group,
    Age = age,
    Gender = gender,
    Weight = weight,
    Temp = vs_24h_temp_first,
    heart_rate = vs_24h_heart_rate_first,
    resp_rate  = vs_24h_resp_rate_first,
    WBC = lab_fst_24h_wbc_first,
    hemoglobin = lab_fst_24h_hemoglobin_first,
    HCT = lab_fst_24h_hct_first,
    PLT = lab_fst_24h_platelet_first,
    BUN = lab_fst_24h_bun_first,
    Cr  = lab_fst_24h_creatinine_first,
    `K+`  = lab_fst_24h_k_first,
    `Na+` = lab_fst_24h_na_first,
    `Cl-` = lab_fst_24h_cl_first,
    PT  = lab_fst_24h_pt_first,
    lactate = lab_fst_24h_lactate_first,
    Hypertension = icd_hypertension,
    Diabetese    = icd_diabetes,
    Heart_Failure = icd_hf,
    COPD = icd_copd,
    Stroke = icd_stroke,
    Renal  = icd_renal,
    liver_cirrhosis = icd_lc,
    CAD = icd_cad,
    Mechanical_ventilation_use    = itvtn_24h_vent_tag,
    Renal_Replacement_Therapy_use = itvtn_24h_rrt_tag,
    Sedative_use                  = drug_24h_sedative_tag,
    NE_equiv_AUC  = ne_equiv_auc_0_24h,
    NE_equiv_Mean = ne_equiv_mean_0_24h,
    NE_equiv_Max  = ne_equiv_max_0_24h
  )
exclude_vars <- intersect(c("stay_id", "Group", "time_lm_days", "event_lm"), names(df))
covars <- setdiff(names(df), exclude_vars)
cat_candidates <- c(
  "Gender",
  "Mechanical_ventilation_use",
  "Renal_Replacement_Therapy_use",
  "Sedative_use",
  "Hypertension","Diabetese","Heart_Failure","COPD","Stroke","Renal","liver_cirrhosis","CAD"
)
cat_covar <- intersect(cat_candidates, covars)
num_covar <- setdiff(covars, cat_covar)

cat("\n[Covars] total =", length(covars),
    " | cat =", length(cat_covar),
    " | num =", length(num_covar), "\n")
ck_num_aligned(df, covars, num_covar)
ck_factor_aligned(df, covars, cat_covar)
vars <- covars
if (!dir.exists("./data")) dir.create("./data", recursive = TRUE)

write.csv(df, file = "./data/map.csv", row.names = FALSE)

save(
  df,
  vars, covars,
  cat_covar, num_covar,
  col_name_gp, gp_name,
  header_name, header_loc,
  file = "./data/df.Rdata"
)

cat("\nDone. Output saved to:\n  ./data/map.csv\n  ./data/df.Rdata\n")
format_num_cell <- function(x, digits = 2) {
  if (all(is.na(x))) return("")
  sprintf(
    paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
    stats::median(x, na.rm = TRUE),
    stats::quantile(x, 0.25, na.rm = TRUE),
    stats::quantile(x, 0.75, na.rm = TRUE)
  )
}

format_cat_cell <- function(x, level) {
  den <- sum(!is.na(x))
  if (den == 0) return("")
  num <- sum(x == level, na.rm = TRUE)
  sprintf("%d (%.1f%%)", num, 100 * num / den)
}

smd_numeric_vs_ref <- function(x, group, ref = "Q2") {
  out <- vapply(setdiff(levels(group), ref), function(gp) {
    x1 <- x[group == gp]
    x0 <- x[group == ref]
    s <- sqrt((stats::var(x1, na.rm = TRUE) + stats::var(x0, na.rm = TRUE)) / 2)
    if (!is.finite(s) || s == 0) return(NA_real_)
    abs((mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) / s)
  }, numeric(1))
  if (all(is.na(out))) NA_real_ else max(out, na.rm = TRUE)
}

smd_categorical_vs_ref <- function(x, group, ref = "Q2") {
  levs <- levels(factor(x))
  out <- c()
  for (gp in setdiff(levels(group), ref)) {
    for (lv in levs) {
      p1 <- mean(x[group == gp] == lv, na.rm = TRUE)
      p0 <- mean(x[group == ref] == lv, na.rm = TRUE)
      denom <- sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
      if (is.finite(denom) && denom > 0) out <- c(out, abs((p1 - p0) / denom))
    }
  }
  if (length(out) == 0 || all(is.na(out))) NA_real_ else max(out, na.rm = TRUE)
}

digits_for_var <- function(v) {
  if (v %in% c("Age", "charlson", "sofa", "sapsii")) return(0)
  if (v %in% c("map24_mean", "NE_equiv_Mean")) return(3)
  2
}

summary_rows <- list()
for (v in vars) {
  if (!v %in% names(df)) next
  if (v %in% cat_covar) {
    x <- factor(df[[v]])
    for (lv in levels(x)) {
      row <- data.frame(
        Variable = v,
        Level = as.character(lv),
        Overall = format_cat_cell(x, lv),
        Q1 = format_cat_cell(x[df$Group == "Q1"], lv),
        Q2 = format_cat_cell(x[df$Group == "Q2"], lv),
        Q3 = format_cat_cell(x[df$Group == "Q3"], lv),
        Q4 = format_cat_cell(x[df$Group == "Q4"], lv),
        SMD_max_vs_Q2 = smd_categorical_vs_ref(x, df$Group),
        stringsAsFactors = FALSE
      )
      summary_rows[[length(summary_rows) + 1]] <- row
    }
  } else {
    x <- suppressWarnings(as.numeric(df[[v]]))
    digs <- digits_for_var(v)
    row <- data.frame(
      Variable = v,
      Level = "",
      Overall = format_num_cell(x, digs),
      Q1 = format_num_cell(x[df$Group == "Q1"], digs),
      Q2 = format_num_cell(x[df$Group == "Q2"], digs),
      Q3 = format_num_cell(x[df$Group == "Q3"], digs),
      Q4 = format_num_cell(x[df$Group == "Q4"], digs),
      SMD_max_vs_Q2 = smd_numeric_vs_ref(x, df$Group),
      stringsAsFactors = FALSE
    )
    summary_rows[[length(summary_rows) + 1]] <- row
  }
}

baseline_summary <- dplyr::bind_rows(summary_rows)
write.csv(baseline_summary, "./data/table1_summary.csv", row.names = FALSE, na = "")
write.csv(baseline_summary, "./data/tbl_smd_body.csv", row.names = FALSE, na = "")
write.csv(
  data.frame(
    Variable = c("Overall", "Q1", "Q2", "Q3", "Q4", "SMD_max_vs_Q2"),
    Description = c(
      "All patients",
      "TWA-MAP quartile 1",
      "TWA-MAP quartile 2",
      "TWA-MAP quartile 3",
      "TWA-MAP quartile 4",
      "Maximum absolute standardized mean difference versus Q2"
    ),
    stringsAsFactors = FALSE
  ),
  "./data/tbl_smd_header.csv",
  row.names = FALSE
)
saveRDS(baseline_summary, "./data/tbl_smd.rds")
cat("Saved reproducible baseline summary to ./data/table1_summary.csv\n")

strong_table_export <- identical(Sys.getenv("MIMICIV_STRONG_TABLE_EXPORT", "0"), "1")
if (strong_table_export) {
tbl_num = 0
cohort_attrib = 'original'
str(df)

table(df$Group)
vars
col_name_gp
vars_to_remove <- c()
for (v in vars) {
  if (all(is.na(df[[v]]))) {
    vars_to_remove <- c(vars_to_remove, v)
    cat("Removed all-NA variable:", v, "\n")
  }
}
if (length(vars_to_remove) > 0) {
  vars <- setdiff(vars, vars_to_remove)
  covars <- setdiff(covars, vars_to_remove)
  num_covar <- setdiff(num_covar, vars_to_remove)
  cat_covar <- setdiff(cat_covar, vars_to_remove)
}
cat("Final number of variables used for Table 1:", length(vars), "\n")
add_tbl_obj(tbl_num)
ok_tbl1(vars, col_name_gp, df, overall=T, method = 'all') 

ok_tbl1(vars, col_name_gp, df, overall=T,smd=F, method = 'auto') 
ok_icm_tbl1(
  vars, col_name_gp, df,smd = T,overall = T,
  var_header_name = header_name,
  var_header_loc = header_loc,method = 'auto'
) %>%
  assign(tbl_name, ., envir = .GlobalEnv)
get(tbl_name)

paste0('Supplementary Table ', update_obj(tbl_num), '. Basic demographic characteristics of the ', cohort_attrib, ' cohort') %>% 
  assign(tbl_cap_name, ., envir = .GlobalEnv)

save_tbl_obj(tbl_name, tbl_cap_name)
docx_status <- tryCatch(
  {
    write_tbl_to_docx(
      write_path = "./data/tbl.docx",
      tbl_name = tbl_name,
      cap_name = tbl_cap_name
    )
    TRUE
  },
  error = function(e) {
    warning("DOCX export skipped: ", conditionMessage(e), call. = FALSE)
    FALSE
  }
)
cat("DOCX export status:", docx_status, "\n")
tbl_num <- 0
add_tbl_obj(tbl_num)

tbl_obj <- ok_icm_tbl1(
  vars,
  col_name_gp,
  df,
  smd = TRUE,
  overall = TRUE,
  var_header_name = header_name,
  var_header_loc = header_loc,
  method = "auto"
)

assign(tbl_name, tbl_obj, envir = .GlobalEnv)
saveRDS(tbl_obj, "./data/tbl_smd.rds")
write.csv(tbl_obj$body$dataset, "./data/tbl_smd_body.csv", row.names = FALSE, na = "")
write.csv(tbl_obj$header$dataset, "./data/tbl_smd_header.csv", row.names = FALSE, na = "")

print(class(tbl_obj))
print(names(tbl_obj))
print(tbl_obj)

cat("Saved SMD Table 1 object to ./data/tbl_smd.rds\n")
cat("Saved SMD Table 1 body to ./data/tbl_smd_body.csv\n")
cat("Saved SMD Table 1 header to ./data/tbl_smd_header.csv\n")
} else {
  cat("Strong formatted table export skipped; set MIMICIV_STRONG_TABLE_EXPORT=1 to enable it.\n")
}

try(DBI::dbDisconnect(con), silent = TRUE)
