# TWA-MAP Sepsis MIMIC-IV: Public Reproduction Code

This repository contains the English public-release SQL and R workflow used to reproduce the downstream extraction, exposure construction, statistical analyses, tables, and figures for the TWA-MAP sepsis study.

## Repository structure

```text
TWA_MAP_Sepsis_MIMICIV/
|-- README.md
|-- README_for_R.md
|-- LICENSE
|-- .gitignore
|-- docs/
|   |-- README_R_Execution_EN.md
|   `-- 00_R_Public_File_Guide_EN.md
|-- scripts/
|   |-- 00_environment_and_directory_check.R
|   |-- 01_main_analysis_generate_cache_and_primary_results.R
|   |-- 02_table1_and_baseline_tables.R
|   |-- 03_major_revision_supplementary_analyses.R
|   |-- 04_q4_deep_phenotype_and_figure1.R
|   |-- 05_observed_mean_no_MI_sensitivity.R
|   |-- 05_plot_preview_in_rstudio.R
|   `-- 06_run_all_in_sequence.R
`-- sql/
    |-- README_SQL_Execution_EN.md
    |-- 01_MAP_trajectory_extraction_EN.sql
    `-- 02_vasopressor_extraction_EN.sql
```

## Scope

This release contains the public SQL and R scripts required to reproduce the final downstream analysis chain from the derived Sepsis-3 source cohort. Historical scripts, exploratory scripts, backup scripts, patient-level data, and generated output files are intentionally excluded.

The SQL workflow starts from an existing Sepsis-3 base cohort and builds the downstream MAP and vasopressor exposure tables used by the R analyses. The upstream derived Sepsis-3 source cohort was constructed using the operational criteria described in the manuscript Methods.

## Requirements

- R with the packages used in the scripts.
- Local PostgreSQL DSN: `mimic4_v31`.
- Authorized access to MIMIC-IV v3.1.
- Either this standalone repository layout or the original archive-style directory structure expected by the supplied scripts.

Database connection settings can be supplied through environment variables:

- `MIMICIV_DSN` (default: `mimic4_v31`)
- `MIMICIV_DATABASE` (default: same as `MIMICIV_DSN`)
- `MIMICIV_DB_SERVER` (default: `localhost`)
- `MIMICIV_DB_PORT` (default: `5432`)
- `MIMICIV_DB_UID` and `MIMICIV_DB_PWD` (optional; omit them when the DSN handles authentication)

## Quick start

Run the SQL scripts first:

1. `sql/01_MAP_trajectory_extraction_EN.sql`
2. `sql/02_vasopressor_extraction_EN.sql`

Then run the R scripts in this order:

1. `scripts/01_main_analysis_generate_cache_and_primary_results.R`
2. `scripts/02_table1_and_baseline_tables.R`
3. `scripts/03_major_revision_supplementary_analyses.R`
4. `scripts/04_q4_deep_phenotype_and_figure1.R`
5. `scripts/05_observed_mean_no_MI_sensitivity.R`

## Data access

This release does not include patient-level data. Access to MIMIC-IV requires completion of PhysioNet credentialing and acceptance of the relevant data use agreement.

## Notes

This repository provides code only. It does not redistribute MIMIC-IV data, patient-level derived tables, local database credentials, or generated analysis outputs.
