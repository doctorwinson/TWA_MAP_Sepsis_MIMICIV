# Supplementary clinical context and follow-up analyses

These scripts accompany the current sepsis landmark manuscript. They extend the
retained F55 analysis (7,289 patients, 948 deaths and 40 completed primary
datasets) without changing cohort eligibility, exposure imputation or the
primary association model. The source scripts are unchanged from those used
to generate the retained supplementary outputs.

## Prerequisites

- Authorized MIMIC-IV v3.1 access and the dependencies in the main
  [workflow README](../../README.md).
- A completed primary workflow with `restricted_cache/eligible_primary.rds`,
  `restricted_cache/analysis_m40.rds`, and `restricted_cache/analysis_specs.rds`.
- R and packages `data.table`, `survival`, and `mice` for `context_revision.R`.
  `followup_summary.R` uses base R and the cached cohort/completed datasets.

Patient-level source data, these RDS files, and completed patient-level
imputation datasets are not distributed. They must be generated locally by
credentialed researchers under PhysioNet terms. Keep credentials outside the
repository. Do not commit restricted caches or database exports.

## Run after the primary workflow

From the `validation_f55` directory, run:

```sh
Rscript extensions/clinical_context/context_revision.R . 02_Results/context_extension
Rscript extensions/clinical_context/followup_summary.R restricted_cache 02_Results/context_extension
```

The first script takes the **analysis folder**, whose `restricted_cache`
subfolder holds its inputs, and creates the output directory. The second
script takes the **restricted-cache folder itself** and requires the output
directory to exist. Do not interchange these arguments. When using other
paths, quote paths containing spaces.

The primary extraction/imputation instructions remain in the main workflow
README; these two commands are not a replacement for that workflow.

## Manuscript output map

| Script/output | Use |
| --- | --- |
| `context_revision.R`: `onset_recording_summary.csv` | Tables S14 and S15: derived onset timing and observed-hour completeness groups |
| `context_revision.R`: `death_location.csv` | Table S16: recorded death timing relative to the index ICU stay, not cause of death |
| `context_revision.R`: `table1_context_rows.csv` | Main Table 1: ICU, admission and advanced-monitoring rows |
| `context_revision.R`: `piecewise_comparison.csv` | Table S10: pooled early/late Q4 contrast from the already specified model |
| `context_revision.R`: `piecewise_estimates_by_imputation.csv` | Aggregate coefficient and variance estimates for the 40 fits |
| `context_revision.R`: `R_sessionInfo.txt` | Runtime/package record; inspect local paths before sharing |
| `followup_summary.R`: `followup_summary.csv` | Results and Table S3: post-landmark follow-up summary |
| `followup_summary.R`: `followup_verification.txt` | Agreement of cohort/event/follow-up fields across completed datasets |

`context_revision.R` refits the previously specified piecewise model to
estimate the direct late-versus-early contrast. Neither script reruns
imputation or searches for a new exposure model. Descriptive quartiles use
the first primary completed dataset, as documented in the manuscript.

## Checks and scope

The scripts contain cohort-size, death-count, identity and expected-value
assertions. Both scripts were syntax parsed before this public update and
compared byte-for-byte with the retained local source. The original scripts
had been run for the retained manuscript outputs; the full extraction,
imputation and analysis pipeline was **not rerun for this code-publication
update**. A fresh one-command end-to-end reproduction has not been certified.

The existing standard MIMIC-derived concepts remain dependencies. Their full
raw-data rebuild was not performed. This code publication does not convert
internal subgroup/period analyses into external validation or establish a
causal mechanism. The complete-recording and later-period weak estimates
remain part of the reported findings.
