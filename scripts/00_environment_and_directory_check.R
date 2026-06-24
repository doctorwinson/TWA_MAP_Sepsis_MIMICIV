# ============================================
# 00. Environment and directory check
# Purpose:
# 1. Detect the current script location and locate the repository/archive root
# 2. Set the working directory to the repository/archive root
# 3. Check required R packages, project structure, and the PostgreSQL DSN
# 4. Create the data/ and outputs_abp_map/ folders if needed
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
archive_root <- find_archive_root(start_dir)
setwd(archive_root)

message("Working directory set to: ", archive_root)

archive_required_dirs <- c(
  "01_SQL\u6570\u636e\u63d0\u53d6",
  "02_R\u7edf\u8ba1\u590d\u73b0",
  "03_\u56fe\u8868",
  "04_\u6587\u7ae0",
  "05_\u9644\u4ef6"
)

repo_required_dirs <- c("scripts", "sql")

is_archive_root <- all(dir.exists(archive_required_dirs))
is_repo_root <- all(dir.exists(repo_required_dirs))
if (!is_archive_root && !is_repo_root) {
  missing_archive_dirs <- archive_required_dirs[!dir.exists(archive_required_dirs)]
  missing_repo_dirs <- repo_required_dirs[!dir.exists(repo_required_dirs)]
  stop(
    "The current working directory is not a valid repository/archive root. ",
    "Missing archive directories: ", paste(missing_archive_dirs, collapse = ", "),
    "; missing repository directories: ", paste(missing_repo_dirs, collapse = ", ")
  )
}

required_pkgs <- c(
  "DBI", "odbc", "data.table", "dplyr", "tidyr",
  "survival", "mice", "splines", "ggplot2",
  "strong", "png", "grid"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  message("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
  message("Install the missing packages before running the downstream scripts.")
} else {
  message("Required R packages are available.")
}

dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs_abp_map", showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("outputs_abp_map", "cache"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("outputs_abp_map", "revision_major"), showWarnings = FALSE, recursive = TRUE)

if (requireNamespace("odbc", quietly = TRUE)) {
  dsn_tbl <- tryCatch(odbc::odbcListDataSources(), error = function(e) NULL)
  if (!is.null(dsn_tbl)) {
    message("Visible ODBC data sources:")
    print(dsn_tbl)
    if (!"mimic4_v31" %in% unlist(dsn_tbl, use.names = FALSE)) {
      message("Note: mimic4_v31 was not found in the current DSN list. Confirm the database configuration before proceeding.")
    }
  }
}

message("Environment check completed. Recommended next step: 01_main_analysis_generate_cache_and_primary_results.R")
