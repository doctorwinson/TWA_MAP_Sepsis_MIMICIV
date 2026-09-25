# F54: corrected landmark, clinical-context explanation and internal stability

This code-only release supersedes prior primary results. It is not an independent external validation and does not claim to identify a causal blood-pressure target or biological mediator.

## Input contract and access

Credentialed access to MIMIC-IV v3.1 is required: https://physionet.org/content/mimiciv/3.1/ (DOI 10.13026/kpb9-mt58).
Patient-level inputs, completed datasets, model objects, logs and outputs must remain in the approved local environment. None are distributed here.

The workflow consumes the restricted `mimic_enriched_preimputation.rds` produced by the preceding F53 extraction/enrichment code, containing 7,770 original source stays, 24 hourly MAP columns, source outcomes/covariates and context enrichments. Set `F54_SOURCE_RDS` explicitly; the fallback path is specific to the author's local archive. Inspect `00_audit_source.R` for assertions and required columns.

The ODBC DSN `mimic4_v31` must expose credentialed MIMIC-IV v3.1 tables, the standard `mimiciv_derived` schema, and the historical local `bdmcc` tables referenced by the audit scripts. Configure credentials outside scripts. SQL is read-only. This code does **not** reconstruct every historical `bdmcc_population` field from raw MIMIC tables. Weight, MV, sedative/analgesic and some coded-condition/intervention definitions remain upstream dependencies requiring their original generation scripts or verified source documentation. Standard first-day weight has a paired-sample sensitivity analysis.

Sepsis flags/timestamps, hourly MAP, recorded death events, age-inclusive Charlson, SOFA and active RRT have independent comparison scripts. Agreement of those fields does not validate the remaining inherited variables.

## Execution

Tested stage-by-stage with Windows, R 4.5.2, data.table, mice, survival, splines, DBI, odbc, ggplot2, ragg and svglite. Runtime package versions are recorded by the scripts in local QC files. Set the working directory to `validation_f54` before running:

```r
Sys.setenv(F54_SOURCE_RDS = "D:/credentialed-local-directory/mimic_enriched_preimputation.rds")
source("01_Code/run_all.R")
```

The runner invokes Rscript child processes, regenerates local outputs and takes substantial time. Its constituent stages were executed and checked; a fresh, one-command full rerun was not separately repeated for this release. It is Windows-oriented (`Rscript.exe`). Never commit the generated restricted caches or output directories.

## Analysis specification

- Sepsis onset, defined as the later qualifying SOFA and suspected-infection time, must be at or before ICU hour24; patients must remain alive and in ICU and have at least18 observed invasive hourly MAP values.
- Five missing source weights are excluded **before** imputation: 7,006 otherwise eligible, 7,001 analyzed, 917 deaths. No patients are removed after imputation.
- Primary and enriched imputation each use m=40, maxit=30, seed42; sex is categorical. Ten-iteration development results are not canonical. The final run was extended from10 to30 after diagnostics, before manuscript estimates were finalized.
- Fixed right-closed quartile boundaries come from the first completed primary dataset. The code retains full precision across imputations and sensitivity analyses.
- Cox estimates and spline contrasts pool all40 imputations. PH diagnostics are separate across imputations, not averaged into a purported pooled P value. Piecewise sensitivities allow all primary coefficients to vary.
- Clinical-context adjustment, SOFA components, leave-one-variable-out, infection timing, observation-density, enriched MI, weight-source and internal time-period analyses are exploratory, not causal mediation.
- Internal temporal comparisons use observed-hour exposures and date-shift-aware possible calendar intervals. They do not establish external generalisability.

Approximate Q4-versus-Q2 estimates: primary2.07, phenotype1.48, extended1.49. Complete24-hour observations1.08 and later-period extended1.07 were near the null. Weak findings must not be omitted when reporting the stronger full-cohort association.

## Script map

`run_all.R` orders the workflow; `00` audits the source; `01` imputes; `02-03` audit cohort/timing; `04-06` model, diagnose and summarize; `07` exports six statistical figures; `08` fits time-varying covariates; `09-12` audit inherited fields and weight/outcome sensitivities; `13` runs downstream stages; `14` verifies statistics; `15` independently reconstructs hourly MAP from raw charted values.

The main cohort flowchart was generated separately using the global Scientific Cohort Flowchart bundled script. Figure1 geometry and document-generation tooling are not reimplemented in this code release. Do not use a historical Figure1 with F54 cohort counts.

## Release checks

The local canonical run passed241 numerical assertions. All released R files parse. This does not certify missing-at-random, prove biological mechanisms, or resolve undocumented upstream definitions. The accompanying manuscript is a journal-neutral transfer draft pending target-journal and author/source confirmations.
