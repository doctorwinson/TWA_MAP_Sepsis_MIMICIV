# Public R File Guide

This document describes the public English R files included in `public_release_en`.

## 1. `scripts/00_environment_and_directory_check.R`

Checks the archive structure, required R packages, output directories, and the visibility of the PostgreSQL DSN used by the project.

## 2. `scripts/01_main_analysis_generate_cache_and_primary_results.R`

Runs the main 24-hour landmark analysis, performs multiple imputation for hourly MAP values, fits the primary and sensitivity Cox models, and generates the main cache and manuscript figure outputs.

## 3. `scripts/02_table1_and_baseline_tables.R`

Builds the baseline dataset from the main analysis cache and generates the manuscript Table 1 dataset together with the full baseline table including standardized mean differences.

## 4. `scripts/03_major_revision_supplementary_analyses.R`

Runs the additional analyses used in the revision package, including quartile mortality summaries, sequential models, subgroup and interaction analyses, missing-data summaries, and proportional-hazards checks.

## 5. `scripts/04_q4_deep_phenotype_and_figure1.R`

Runs the Q4 deep-phenotype analysis, the phenotype-augmented model, and the final Figure 1 flow diagram.

## 6. `scripts/05_plot_preview_in_rstudio.R`

Reads the final PNG figures already stored in the archive and displays them sequentially in the RStudio `Plots` pane.

## 7. `scripts/06_run_all_in_sequence.R`

Runs the core reproduction chain in sequence by calling files `00` through `04`.

## Notes

- These files are the public English release intended for direct reproduction.
- The original archive may still contain additional internal, historical, or backup files outside this folder.
- The analysis assumes that the required SQL outputs have already been created in the database.
