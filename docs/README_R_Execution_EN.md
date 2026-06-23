# R Execution Guide (Public English Release)

This directory contains the English public-release version of the R workflow required to reproduce the main analyses, tables, and figures.

## Included files

- `scripts/00_environment_and_directory_check.R`
- `scripts/01_main_analysis_generate_cache_and_primary_results.R`
- `scripts/02_table1_and_baseline_tables.R`
- `scripts/03_major_revision_supplementary_analyses.R`
- `scripts/04_q4_deep_phenotype_and_figure1.R`
- `scripts/05_plot_preview_in_rstudio.R`
- `scripts/06_run_all_in_sequence.R`
- `docs/00_R_Public_File_Guide_EN.md`

## Recommended execution order

1. `scripts/00_environment_and_directory_check.R`
2. `scripts/01_main_analysis_generate_cache_and_primary_results.R`
3. `scripts/02_table1_and_baseline_tables.R`
4. `scripts/03_major_revision_supplementary_analyses.R`
5. `scripts/04_q4_deep_phenotype_and_figure1.R`
6. `scripts/05_plot_preview_in_rstudio.R`

To run the core workflow in sequence, use:

- `scripts/06_run_all_in_sequence.R`

## Scope

This public English set includes only the files required for direct reproduction of the final analysis chain. Historical, backup, and exploratory scripts are intentionally excluded from this folder.

## Runtime assumptions

- Run the scripts from RStudio or R in either this standalone repository or the existing archive structure.
- The scripts automatically locate the repository/archive root and set the working directory.
- The local PostgreSQL DSN expected by the scripts is `mimic4_v31` by default.
- Optional connection details can be supplied through `MIMICIV_DSN`, `MIMICIV_DATABASE`, `MIMICIV_DB_SERVER`, `MIMICIV_DB_PORT`, `MIMICIV_DB_UID`, and `MIMICIV_DB_PWD`.
- The repository/archive should already contain the `data/` and `outputs_abp_map/` folders or allow the scripts to create them.
