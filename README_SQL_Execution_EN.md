# SQL Execution Guide

This folder contains the SQL workflow used to build the main analysis tables for the invasive ABP-MAP study.

## Files

- `01_MAP_trajectory_extraction_EN.sql`
- `02_vasopressor_extraction_EN.sql`


## Recommended Execution Order

1. `01_MAP_trajectory_extraction_EN.sql`
2. `02_vasopressor_extraction_EN.sql`

## What Each Script Does

### 1. `01_MAP_trajectory_extraction_EN.sql`

This script builds the main 24-hour landmark cohort and extracts invasive arterial MAP measurements during the first 24 ICU hours.

It creates:

- `bdmcc.cohort_design_abp`
- `bdmcc.map_hourly_final_abp`
- `bdmcc.base_landmark_surv_abp`

Main features:

- starts from `bdmcc.bdmcc_population`
- applies the sepsis cohort restrictions already encoded in that table
- defines a 24-hour landmark after ICU admission
- flags pre-landmark death records
- extracts invasive ABP-MAP only (`itemid = 220052`)
- aggregates MAP values by ICU hour
- applies the hourly coverage threshold for inclusion

## 2. `02_vasopressor_extraction_EN.sql`

This script extracts vasopressor exposure during the first 24 ICU hours and converts drug-specific exposure to norepinephrine-equivalent values.

It creates:

- `bdmcc.vaso_ne_equiv_0_24h`

Main features:

- reads infusion records from `mimiciv_icu.inputevents`
- clips exposure to the first 24 ICU hours
- standardizes supported drug units
- calculates duration and AUC for each infusion segment
- derives 24-hour norepinephrine-equivalent mean, max, and AUC

## Required Upstream Tables

Before running these scripts, the following tables must already be available:

- `bdmcc.bdmcc_population`
- `mimiciv_hosp.admissions`
- `mimiciv_icu.icustays`
- `mimiciv_icu.chartevents`
- `mimiciv_icu.inputevents`
- `mimiciv_icu.d_items`

## Important Scope Note

This archive does **not** include the SQL used to build `bdmcc.bdmcc_population` from raw MIMIC-IV tables.

The current SQL workflow therefore starts from an existing Sepsis-3 base cohort and covers the downstream extraction steps used for the present study.


These R scripts perform:

- hourly MAP processing
- multiple imputation
- survival modeling
- sensitivity analyses
- figure and table generation

## Quality Control

Both SQL scripts include post-run QC queries at the end. These checks should be reviewed after execution to confirm that:

- output tables were created successfully
- cohort counts are plausible
- no ineligible records remain in the final landmark cohort
- norepinephrine-equivalent exposure values are populated as expected

## Suggested Practice

- Run the scripts in a PostgreSQL environment connected to MIMIC-IV v3.1.
- Keep the execution log for record keeping.
- Export the final tables to CSV after successful execution if a local archive is required.
