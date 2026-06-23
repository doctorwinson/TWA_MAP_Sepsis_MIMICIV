# ============================================
# 00. Environment and directory check
# Purpose:
# 1. Detect the current script location and locate the archive root
# 2. Set the working directory to the archive root
# 3. Check required R packages, archive structure, and the PostgreSQL DSN
# 4. Create the data/ and outputs_abp_map/ folders if needed
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
    if (all(dir.exists(file.path(current_dir, required_dirs)))) {
      return(current_dir)
    }
    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      stop("Unable to locate the archive root automatically. Run this script from the archive root or one of its subdirectories.")
    }
    current_dir <- parent_dir
  }
}

this_script <- locate_current_script()
start_dir <- if (!is.na(this_script) && nzchar(this_script)) dirname(this_script) else getwd()
archive_root <- find_archive_root(start_dir)
setwd(archive_root)

message("Working directory set to: ", archive_root)

required_dirs <- c(
  "01_SQL数据提取",
  "02_R统计复现",
  "03_图表",
  "04_文章",
  "05_附件"
)

missing_dirs <- required_dirs[!dir.exists(required_dirs)]
if (length(missing_dirs) > 0) {
  stop("The current working directory is not a valid archive root. Missing directories: ", paste(missing_dirs, collapse = ", "))
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
