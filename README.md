# TWA-MAP Sepsis MIMIC-IV: Public Reproduction Code

## Current reconstructed workflow: F55 (2026-09-25)

Use [validation_f55](validation_f55/README.md) for the explicit study-specific upstream reconstruction and current analyses: 7,289 patients,948 deaths; two m=40/30-iteration imputations. Historical bdmcc cohort fields are no longer required, while documented standard MIMIC-derived concepts remain dependencies. Weak complete-observation and later-period results are retained. No patient-level data are released.

## Superseded corrected-cohort workflow: F54

Use [validation_f54](validation_f54/README.md) for the corrected hour-24 sepsis-landmark analysis. It supersedes the legacy primary results: 7,001 patients, 917 deaths, two imputation specifications with 40 datasets and 30 iterations each. The public release contains code only, not patient-level data. Upstream source dependencies remain explicitly documented; this is not a raw-MIMIC-only reconstruction of every inherited variable. Internal calendar-period comparisons are not external validation, and complete-observation/later-period results were weak.

The original Scientific Reports submission is closed. Older scripts below and `validation_f53` are retained for provenance and must not be mistaken for the current F54 results.

## Legacy workflow

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

The legacy scripts reproduce the earlier downstream analysis chain from the derived Sepsis-3 source cohort. They are retained as historical code; they are not the canonical corrected F54 analysis. Patient-level data and generated outputs are excluded from this repository.

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

## Post hoc validation extension

The separate [`validation_f53`](validation_f53/README.md) package adds clinical-context and SOFA analyses, retrospective internal temporal stability checks, and exploratory eICU in-hospital mortality transportability. It preserves the original workflow and explicitly documents source-cohort dependencies, landmark-onset eligibility findings, missing-data conventions and the difference between hospital and day-30 outcomes. It does not claim definitive external validation or causal pathway identification.

## Notes

This repository provides code only. It does not redistribute MIMIC-IV data, patient-level derived tables, local database credentials, or generated analysis outputs.
