# ============================================
# 06. Run the full public workflow in sequence
# Purpose:
# 1. Locate the archive root automatically
# 2. Run the environment check and the core analysis scripts in order
# 3. Keep the execution chain suitable for direct use in RStudio
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

message("Working directory set to: ", getwd())

run_scripts <- c(
  "00_environment_and_directory_check.R",
  "01_main_analysis_generate_cache_and_primary_results.R",
  "02_table1_and_baseline_tables.R",
  "03_major_revision_supplementary_analyses.R",
  "04_q4_deep_phenotype_and_figure1.R"
)

for (script_name in run_scripts) {
  script_path <- file.path("02_R统计复现", "public_release_en", "scripts", script_name)
  message(">>> Running: ", script_path)
  source(script_path, encoding = "UTF-8", local = new.env(parent = globalenv()))
  message("<<< Completed: ", script_path)
}

message("The public reproduction workflow completed successfully.")
