# Post hoc phenotype explanation and validation extension

This directory adds an exploratory extension without replacing the original study scripts. It is a research analysis package, not an independently validated clinical model.

## Scope

- Reconstruct the complete primary m=40 MICE object from the preserved pre-imputation analysis cache.
- Compare clinical-context adjustment blocks and SOFA components on fixed analytic samples.
- Assess retrospective internal temporal stability using admission-year intervals rather than shifted calendar years.
- Report ICU/service heterogeneity, interaction tests, early-sepsis eligibility sensitivity and complete-observation sensitivity.
- Compare in-hospital mortality associations with a separately analysed eICU cohort. This is not external validation of complete day-30 vital status.
- Save diagnostics and aggregate figure source data locally.

## Required local inputs

The code contains no patient data. Credentialed local access and existing source tables are required.

1. Set `F53_SOURCE_ROOT` to a directory containing `outputs_abp_map/cache/cache_main_ge18h_Model2_SOFA_minimal.rds`, produced by the original workflow. This must preserve the original pre-imputation 7,770-patient cohort. The old first-imputation cache is not treated as a complete m=40 object.
2. Configure the `mimic4_v31` and `eicu` ODBC DSNs using local authentication. New database queries are read-only.
3. MIMIC-IV requires the original `bdmcc` source cohort and the standard hospital, ICU and derived SOFA tables.
4. eICU requires `bdmcc.bdmcc_population`, `bdmcc.sepsis3`, `eicu_icu.patient` and `eicu_icu.vitalperiodic`.

**The local eICU Sepsis-3 source-table creation SQL was not available for verification in this extension.** This package does not reconstruct that definition from raw eICU tables. Its antibiotic/culture/SOFA timing field is a proxy sensitivity, not an independently verified strict Sepsis-3 definition. Pregnancy exclusion and first-ever ICU admission are not fully harmonized. Do not claim end-to-end external reproducibility until those dependencies are resolved.

## Run

From this directory:

```r
Sys.setenv(F53_SOURCE_ROOT = "/path/to/authorized/source")
source("01_Code/run_all.R")
```

Alternatively run `Rscript 01_Code/run_all.R` after setting the environment variable. Dependencies include R, data.table, mice, mitml, survival, DBI, odbc, sandwich, ggplot2 and optional svglite. No automatic package installation or database mutation is included.

Outputs are written to ignored `restricted_cache`, `02_Results` and `04_QC` directories. Restricted caches include imputed patient-level datasets and fitted models; they must not be committed or redistributed through this repository. The plotting script emits PDF, SVG (when available), PNG and TIFF files, with 600 dpi raster exports.

## Interpretation and newly identified limitations

The primary reconstructed estimate and phenotype-adjusted estimate match the submitted version. The source cohort includes patients whose recorded sepsis onset is after the 24-hour landmark; this is explicitly audited, with an early-onset restriction. A redesigned primary cohort would require fresh cohort-specific imputation and synchronized manuscript changes.

Internal periods were already included in discovery, so this is retrospective stability rather than independent held-out validation. Adjusted-coefficient changes are not mediated effects or causal contribution percentages. SOFA domain analyses distinguish source-convention missing-as-zero scoring from complete raw-domain cases. The original character-sex imputation-predictor omission is preserved for reproduction and addressed in phenotype-enriched imputation sensitivities. Diagnostics retain constant/redundant predictor omissions.

All extensions are post hoc. Effect heterogeneity, weak complete-observation results and endpoint differences must be reported, irrespective of statistical significance.

See `analysis_plan.md` and `amendments.md` for the analysis specification and audit-driven changes.
