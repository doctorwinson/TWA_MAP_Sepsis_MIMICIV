# F55: explicit upstream reconstruction and observational stability analyses

This is the current study-specific analysis workflow. It supersedes F54 for the proposed PLOS ONE manuscript. It is not an external validation, causal mediation analysis or an estimate of an optimal MAP target.

## Data access and dependencies

Credentialed MIMIC-IV v3.1 access is required: https://physionet.org/content/mimiciv/3.1/ (DOI: 10.13026/kpb9-mt58). No patient-level data, fitted model objects or imputed datasets are distributed.

The read-only PostgreSQL connection uses ODBC DSN `mimic4_v31` by default; set `MIMIC_DSN` to override. Configure credentials outside the scripts. Required schemas are `mimiciv_hosp`, `mimiciv_icu` and standard `mimiciv_derived` concepts: sepsis3, first_day_weight, charlson, first_day_sofa, first_day_rrt, sapsii and ventilation. Study-specific extraction no longer requires historical `bdmcc` cohort tables.

The standard derived concepts remain dependencies; rebuilding their entire schema from raw data was not performed. `standard_concept_reference/manifest.json` identifies pinned official reference definitions, not the historical build version of the local derived schema. It does not certify that all installed concepts match that commit.

## Execution

Windows R 4.5.2; packages DBI, odbc, data.table, mice, survival, splines, ggplot2, ragg and svglite. Set the working directory to `validation_f55`, then run:

```r
source("01_Code/run_all.R")
```

The runner creates local output/cache directories, extracts the source, runs two m=40/30-iteration imputations and regenerates analyses, figures and independent checks. The constituent stages were run and verified locally; a separate fresh one-command end-to-end reproduction was not repeated. The runner uses Windows `Rscript.exe`. Output and restricted cache folders are ignored by Git. Do not upload them or credentials.

Optional `F55_LEGACY_AUDIT_RDS` supplies an authorized F54 cache for a local aggregate change audit. Its absence does not affect the primary analysis. Do not set `F55_SKIP_EXTRACT=1` for fresh reproduction. Historical `F54_MI_SCOPE`, `F54_MAXIT` and `F54_RESUME_MI` environment-variable names remain for compatibility; the runner sets scope and iterations explicitly and disables resume.

## Reconstruction and analysis

- Eligibility: adult first recorded ICU stay, still in ICU and alive after hour24, retrospectively defined Sepsis-3 onset by24h, no specified malignancy/pregnancy codes, at least18 observed invasive hourly MAPs.
- Final flow: 7,305 with sufficient observation; 16 missing first-day weights excluded before MI; 7,289 analyzed and 948 deaths. No post-MI exclusions.
- MAP: raw chartevents220052, 20-200mmHg, hourly means in [0,24h); all24 hours after predictive mean matching. Two specifications use m40, maxit30, seed42.
- Weight: standard first-day mean weight; invasive ventilation: overlapping InvasiveVent interval, not tracheostomy alone; sedation/analgesia: explicit positive administrations in inputevents; rates and clipped infusion intervals audited.
- Hospitalization ICD codes do not establish pre-existing disease. CABG/PCI are hospitalization procedure codes, not proven early interventions. Legacy column names containing `24h_cabg` or `24h_pci` are compatibility names only. Standard SOFA/RRT windows can extend6h before ICU entry.
- Malignancy includes ICD9 140-208,2090-2093,2097 and ICD10 C. Benign2094-2096 are not excluded as malignancy.
- Death follow-up runs from ICU hour24 to admission day30 (maximum29days); hospital death timestamps take precedence over date-only linked death records.
- Fixed quartile cutpoints come from the first primary completed dataset; all40 regressions pool under Rubin's rules. Descriptive tables use the first dataset.
- Complete-observation, observed-hour, stricter infection-pair, alternative definitions, SOFA-component, clinical-context and temporal analyses are exploratory. No models were added to obtain significance.

Final Q4-versus-Q2 HRs (95% CI): primary2.10(1.72-2.56), phenotype1.48(1.21-1.81), extended1.47(1.20-1.81), complete24h1.12(0.77-1.63), later-period extended1.10(0.81-1.49). The boundary-period extended model (274 patients,27 deaths) was not reliably estimable; its failure is explicitly retained in output. These results must not be mixed with legacy versions.

## Verification and script map

`19` reconstructs; `01` imputes; `04` models; `05` splines/diagnostics; `06` enriched MI/tables; `07` figures; `08` time-varying coefficients; `12` paired definition sensitivities; `13` orders downstream stages; `22` independently reaggregates raw hourly MAP/death; `25` rechecks weight/RRT/SOFA; `14` verifies results. The local final run passed238 numerical assertions; raw MAP/outcome discrepancies were zero.

All released R scripts are syntax checked. These checks do not establish missing-at-random or clinical causality. Figure1 uses the separate global Scientific Cohort Flowchart bundled renderer with the included count configuration; geometry and private manuscript tooling are not reimplemented here.
