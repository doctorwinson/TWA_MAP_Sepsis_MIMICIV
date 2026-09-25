# F53 phenotype explanation and validation analysis plan

Locked locally on 2026-09-25 before viewing new association estimates. This is a post-rejection, post hoc extension, not prospective preregistration.

## Objectives and boundaries

1. Reproduce the submitted F52 m=40 primary and phenotype-adjusted analyses from the preserved pre-imputation source, saving the complete imputation object and diagnostics.
2. Assess within-centre temporal stability and clinical-service heterogeneity of the achieved-MAP association.
3. Quantify changes in the Q4 versus Q2 coefficient after adjustment for clinical-context blocks and individual SOFA domains. These are descriptive coefficient changes, not mediated proportions or causal pathways. Correlated blocks, non-collapsibility, measurement overlap and confounding by indication preclude additive attribution.
4. Examine eICU as a separate multicentre transportability cohort. Its routinely recorded outcome is hospital discharge status, not complete day-30 vital status. Use a harmonized in-hospital endpoint in both databases for this secondary analysis. Never relabel this as external validation of day-30 mortality.

## Frozen MIMIC-IV definitions

- Adult Sepsis-3 cohort and 24-hour landmark definitions remain those of the archived F52 cohort. Preserve its original selection for reproduction; separately report sepsis onset after the landmark and evaluate restriction to onset by 24 h if present.
- Primary day-30 endpoint and follow-up use the preserved time_lm_days and event_lm, with an audit against source dates.
- Rebuild PMM MICE m=40, 10 iterations, seed 42 using the original function, adjacent hourly values, primary covariates and outcome predictors. Record all excluded rows and residual missing covariates. Pool model coefficients and variances using Rubin's rules and Barnard-Rubin degrees of freedom.
- Calculate exact quartile thresholds from the first rebuilt imputation as in the source code. Compare with F52 rounded thresholds 71.40, 75.72 and 80.99 mmHg. Fixed thresholds across imputations. External observed-hour analyses use the published rounded thresholds, with infinite outer boundaries, and clearly identify absolute categories rather than local quartiles.
- Primary adjustment: age, sex, weight, Charlson score, total SOFA, MV, RRT and sedative use during 0-24 h.

## Temporal and service analyses

- Derive the possible real admission-year interval by adding year(intime)-anchor_year to both bounds of anchor_year_group. Do not split using de-identified calendar years alone.
- Earlier cohort: latest possible admission year <=2016. Later cohort: earliest possible admission year >=2017. Exclude intervals straddling this boundary from the primary temporal comparison and enumerate them. Use observed-hour mean MAP and frozen F52 categories for the primary temporal stability analysis, avoiding cross-period imputation. Complete 24/24-hour restriction is a sensitivity analysis.
- These periods were already included in the original discovery cohort. Call this retrospective internal temporal stability, not independent held-out validation.
- Clinical-context strata: medical, cardiac, surgical/trauma, neuro ICU/stepdown and other; additionally identify hospital service at ICU entry from the most recent services record at or before intime. ICU-type strata alone are within-centre heterogeneity analyses, not external validation.
- Report events, denominators, HRs and CIs; suppress unstable subgroup coefficients rather than interpreting significance. Interaction tests use common analytic samples and pooled models where applicable. Multiplicity-adjust exploratory interaction p-values using Benjamini-Hochberg.

## Phenotype blocks

- Chronic vascular context: hypertension and prior stroke codes, recognizing coding may not establish onset before admission.
- ICU and admission context: ICU type, admission type, service at ICU entry.
- Cardiovascular procedures/devices: CABG, PCI, IABP; advanced monitoring: PICCO/NICOM.
- Support intensity: norepinephrine-equivalent mean dose (log1p sensitivity), MV, RRT, sedatives. MV/RRT/sedatives are already in the primary model; evaluate leave-one-block-out and base-to-augmented models without double counting them.
- SOFA components: replace the total score with six domain scores; also replace it with non-cardiovascular SOFA and with non-CNS SOFA. Never adjust for total SOFA and all its components in the same model. Compare the source total and domain sum, and keep equal samples for coefficient comparisons.
- Exclusion analyses: remove neuro ICU, surgical/trauma ICU, or cardio-device cases, plus a medical-only analysis, using the same thresholds. Findings remain observational.
- No formal mediation, causal direct effect or mechanistic percentage is claimed because MAP, organ dysfunction and early interventions overlap in time.

## eICU feasibility and secondary transportability

- Audit existing bdmcc derived Sepsis-3 tables and document their provenance. Existing materialized flags without derivation source are provisional. Inspect antibiotic/culture-plus-SOFA timing separately from diagnosis-driven criteria.
- Adults, first ICU stay in a hospitalization (eICU does not provide exact chronology across all hospitalizations), alive and still in ICU at 24 h, known discharge outcome, no malignancy for aligned sensitivity, >=18 observed arterial-MAP hours in [0,1440) minutes. Account for hospital clustering and repeated patients where relevant. Pregnancy alignment requires explicit ascertainment; unavailable ascertainment is a reported limitation, not an assumed negative.
- Exposure: systemicMean from vitalPeriodic, 20-200 mmHg, mean within each hour, equally weighted mean across observed hours. Complete 24-hour restriction. No pooling with MIMIC-IV before estimation.
- Endpoint: in-hospital mortality. Use a prespecified shared adjustment set available in both datasets. Hospital-stratified Cox is not equivalent to a day-30 all-cause endpoint; prefer logistic OR with hospital-clustered robust uncertainty for this secondary comparison.
- No selection of thresholds, subgroups or outcome definitions based on significance. Report incompatible or null results.

## Deliverables

Versioned scripts/SQL, local restricted caches, aggregate CSV tables, publication-quality figures, editable bilingual review report and an English supplementary analysis document. Patient-level datasets and identifiers must not enter the public code release. All analyses remain post hoc. A new target journal and submission-ready main manuscript will follow interpretation of the new evidence.

## Primary sources checked

- MIMIC-IV v3.1: https://physionet.org/content/mimiciv/3.1/
- Admission-year anchoring: https://github.com/MIT-LCP/mimic-iv/blob/master/concepts/demographics/age.sql
- MIMIC-III overlap and CareVue subset: https://physionet.org/content/mimic3-carevue/1.4/
- eICU outcomes: https://eicu.mit.edu/eicutables/patient/
- eICU arterial-pressure provenance: https://eicu.mit.edu/eicutables/vitalperiodic/
- Multiple-imputation pooling: https://amices.org/mice/reference/pool.html
