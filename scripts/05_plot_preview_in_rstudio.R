# ============================================
# 05. Plot preview in the RStudio Plots pane
# Purpose:
# 1. Locate the repository/archive root automatically
# 2. Read the final PNG figures stored in the project
# 3. Display the manuscript and supplementary figures in sequence
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
  archive_dirs <- c("01_SQL数据提取", "02_R统计复现", "03_图表", "04_文章", "05_附件")
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

suppressPackageStartupMessages({
  library(png)
  library(grid)
})

fig_paths <- c(
  "Figure 1" = file.path("03_图表", "01_主文图片", "Figure_1_flow_diagram.png"),
  "Figure 2" = file.path("03_图表", "01_主文图片", "Figure_2_TWA_MAP_spline.png"),
  "Figure 3" = file.path("03_图表", "01_主文图片", "Figure_3_TTRlt65_spline.png"),
  "Figure S1" = file.path("03_图表", "02_补充图片", "Supplementary_Figure_S1_ge20h_spline.png"),
  "Figure S2" = file.path("03_图表", "02_补充图片", "Supplementary_Figure_S2_complete24h_spline.png"),
  "Figure S3" = file.path("03_图表", "02_补充图片", "Supplementary_Figure_S3_SAPSII_spline.png"),
  "Figure S4" = file.path("03_图表", "02_补充图片", "Supplementary_Figure_S4_NE_adjusted_spline.png")
)

show_png <- function(path, title_text) {
  if (!file.exists(path)) {
    warning("File not found: ", path)
    return(invisible(NULL))
  }
  img <- png::readPNG(path)
  grid::grid.newpage()
  grid::grid.raster(img)
  grid::grid.text(
    label = title_text,
    x = 0.02, y = 0.98,
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
}

message("Working directory set to: ", getwd())
message("Displaying figures in the RStudio Plots pane. Use the Plots history to review them.")

for (nm in names(fig_paths)) {
  show_png(fig_paths[[nm]], nm)
}
