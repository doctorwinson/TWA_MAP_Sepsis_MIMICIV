/* =====================================================================================
   Purpose:
   - Extract vasopressor exposure during the first 24 ICU hours
   - Derive norepinephrine-equivalent exposure for each stay

   Data source:
   - mimiciv_icu.inputevents

   Output table:
   - bdmcc.vaso_ne_equiv_0_24h

   Derived measures:
   - ne_equiv_mean_0_24h
   - ne_equiv_max_0_24h
   - ne_equiv_auc_0_24h

   Unit notes:
   - AUC unit: (ug/kg/min) * hour
   - Conversion factors:
     norepinephrine = 1
     epinephrine    = 1
     dopamine       = 1/100
     phenylephrine  = 1/10
     vasopressin    = 2.5 (when rate is in units/min)
===================================================================================== */

------------------------------------------------------------
-- Step 0: Drop existing output and temporary tables
------------------------------------------------------------
DROP TABLE IF EXISTS bdmcc.vaso_ne_equiv_0_24h;

DROP TABLE IF EXISTS tmp_vaso_raw;
DROP TABLE IF EXISTS tmp_vaso_clip;
DROP TABLE IF EXISTS tmp_vaso_std;
DROP TABLE IF EXISTS tmp_vaso_auc;
DROP TABLE IF EXISTS tmp_ne_equiv;

------------------------------------------------------------
-- Step 1: Get eligible stays and ICU admission time
------------------------------------------------------------
CREATE TEMP TABLE tmp_cohort AS
SELECT
    c.stay_id,
    i.intime
FROM bdmcc.cohort_design_abp c
JOIN mimiciv_icu.icustays i
  ON c.stay_id = i.stay_id
WHERE c.landmark_eligible = 1;

------------------------------------------------------------
-- Step 2: Extract raw vasopressor infusion records overlapping 0-24 hours
------------------------------------------------------------
CREATE TEMP TABLE tmp_vaso_raw AS
SELECT
    ie.stay_id,
    ie.starttime,
    ie.endtime,
    ie.itemid,
    LOWER(di.label) AS label,
    ie.rate,
    LOWER(COALESCE(ie.rateuom,'')) AS rateuom
FROM mimiciv_icu.inputevents ie
JOIN mimiciv_icu.d_items di
  ON ie.itemid = di.itemid
JOIN tmp_cohort tc
  ON ie.stay_id = tc.stay_id
WHERE ie.starttime < tc.intime + INTERVAL '24 hours'
  AND COALESCE(ie.endtime, ie.starttime) > tc.intime
  AND ie.rate IS NOT NULL
  AND ie.rate > 0
  AND (
        di.label ILIKE '%norepinephrine%'
     OR di.label ILIKE '%epinephrine%'
     OR di.label ILIKE '%dopamine%'
     OR di.label ILIKE '%vasopressin%'
     OR di.label ILIKE '%phenylephrine%'
  );

------------------------------------------------------------
-- Step 3: Clip each infusion to the first 24 ICU hours
------------------------------------------------------------
CREATE TEMP TABLE tmp_vaso_clip AS
SELECT
    v.stay_id,
    GREATEST(v.starttime, tc.intime) AS t0,
    LEAST(COALESCE(v.endtime, v.starttime), tc.intime + INTERVAL '24 hours') AS t1,
    v.itemid,
    v.label,
    v.rate,
    v.rateuom
FROM tmp_vaso_raw v
JOIN tmp_cohort tc
  ON v.stay_id = tc.stay_id
WHERE LEAST(COALESCE(v.endtime, v.starttime), tc.intime + INTERVAL '24 hours')
    > GREATEST(v.starttime, tc.intime);

------------------------------------------------------------
-- Step 4: Standardize drug labels and rate units
------------------------------------------------------------
CREATE TEMP TABLE tmp_vaso_std AS
SELECT
    stay_id,
    t0, t1,
    CASE
        WHEN label LIKE '%norepinephrine%' THEN 'norepinephrine'
        WHEN label LIKE '%epinephrine%'    THEN 'epinephrine'
        WHEN label LIKE '%dopamine%'       THEN 'dopamine'
        WHEN label LIKE '%vasopressin%'    THEN 'vasopressin'
        WHEN label LIKE '%phenylephrine%'  THEN 'phenylephrine'
        ELSE 'other'
    END AS drug,

    CASE
        WHEN (label LIKE '%norepinephrine%' OR label LIKE '%epinephrine%'
           OR label LIKE '%dopamine%' OR label LIKE '%phenylephrine%')
         AND (rateuom LIKE '%ug/kg/min%' OR rateuom LIKE '%mcg/kg/min%')
            THEN rate

        WHEN label LIKE '%vasopressin%'
         AND (rateuom LIKE '%units/min%' OR rateuom LIKE '%unit/min%')
            THEN rate

        WHEN label LIKE '%vasopressin%'
         AND (rateuom LIKE '%units/hour%' OR rateuom LIKE '%unit/hour%')
            THEN rate / 60.0

        ELSE NULL
    END AS rate_std
FROM tmp_vaso_clip;

------------------------------------------------------------
-- Step 5: Compute infusion duration and AUC contribution
------------------------------------------------------------
CREATE TEMP TABLE tmp_vaso_auc AS
SELECT
    stay_id,
    drug,
    rate_std,
    EXTRACT(EPOCH FROM (t1 - t0))/3600.0 AS dur_h,
    (rate_std * EXTRACT(EPOCH FROM (t1 - t0))/3600.0) AS auc_part
FROM tmp_vaso_std
WHERE rate_std IS NOT NULL
  AND EXTRACT(EPOCH FROM (t1 - t0)) > 0;

------------------------------------------------------------
-- Step 6: Summarize 0-24h AUC by drug
------------------------------------------------------------
CREATE TEMP TABLE tmp_drug_summary AS
SELECT
    stay_id,
    drug,
    SUM(auc_part) AS auc_0_24h,
    MAX(rate_std) AS rate_max
FROM tmp_vaso_auc
GROUP BY stay_id, drug;

------------------------------------------------------------
-- Step 7: Convert drug-specific exposure to norepinephrine-equivalent values
------------------------------------------------------------
CREATE TEMP TABLE tmp_ne_equiv AS
SELECT
    stay_id,
    SUM(CASE WHEN drug='norepinephrine' THEN auc_0_24h ELSE 0 END) AS ne_auc,
    SUM(CASE WHEN drug='epinephrine'    THEN auc_0_24h ELSE 0 END) AS epi_auc,
    SUM(CASE WHEN drug='dopamine'       THEN auc_0_24h ELSE 0 END) / 100.0 AS dopa_auc_ne,
    SUM(CASE WHEN drug='phenylephrine'  THEN auc_0_24h ELSE 0 END) / 10.0 AS phenyl_auc_ne,
    SUM(CASE WHEN drug='vasopressin'    THEN auc_0_24h ELSE 0 END) * 2.5 AS vaso_auc_ne,

    MAX(CASE WHEN drug='norepinephrine' THEN rate_max ELSE 0 END) AS ne_max,
    MAX(CASE WHEN drug='epinephrine'    THEN rate_max ELSE 0 END) AS epi_max,
    MAX(CASE WHEN drug='dopamine'       THEN rate_max ELSE 0 END) / 100.0 AS dopa_max_ne,
    MAX(CASE WHEN drug='phenylephrine'  THEN rate_max ELSE 0 END) / 10.0 AS phenyl_max_ne,
    MAX(CASE WHEN drug='vasopressin'    THEN rate_max ELSE 0 END) * 2.5 AS vaso_max_ne
FROM tmp_drug_summary
GROUP BY stay_id;

------------------------------------------------------------
-- Step 8: Save the final norepinephrine-equivalent table
------------------------------------------------------------
CREATE TABLE bdmcc.vaso_ne_equiv_0_24h AS
SELECT
    c.stay_id,

    COALESCE(
        ne.ne_auc + ne.epi_auc + ne.dopa_auc_ne + ne.phenyl_auc_ne + ne.vaso_auc_ne,
        0
    ) AS ne_equiv_auc_0_24h,

    COALESCE(
        (ne.ne_auc + ne.epi_auc + ne.dopa_auc_ne + ne.phenyl_auc_ne + ne.vaso_auc_ne) / 24.0,
        0
    ) AS ne_equiv_mean_0_24h,

    GREATEST(
        COALESCE(ne.ne_max, 0),
        COALESCE(ne.epi_max, 0),
        COALESCE(ne.dopa_max_ne, 0),
        COALESCE(ne.phenyl_max_ne, 0),
        COALESCE(ne.vaso_max_ne, 0)
    ) AS ne_equiv_max_0_24h

FROM tmp_cohort c
LEFT JOIN tmp_ne_equiv ne
  ON c.stay_id = ne.stay_id
ORDER BY c.stay_id;

------------------------------------------------------------
-- Step 9: Post-run quality checks
------------------------------------------------------------

-- 1) Confirm that the output table was created
SELECT to_regclass('bdmcc.vaso_ne_equiv_0_24h') AS vaso_table;

-- 2) Count all eligible stays in the vasopressor table
SELECT COUNT(*) AS n_total
FROM bdmcc.vaso_ne_equiv_0_24h;

-- 3) Count stays with and without vasopressor exposure
SELECT
    COUNT(*) AS n_total,
    SUM(CASE WHEN ne_equiv_auc_0_24h > 0 THEN 1 ELSE 0 END) AS n_with_vaso,
    SUM(CASE WHEN ne_equiv_auc_0_24h = 0 THEN 1 ELSE 0 END) AS n_without_vaso
FROM bdmcc.vaso_ne_equiv_0_24h;

-- 4) Descriptive statistics for norepinephrine-equivalent exposure
SELECT
    MIN(ne_equiv_mean_0_24h) AS min_mean,
    PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY ne_equiv_mean_0_24h) AS p25_mean,
    PERCENTILE_CONT(0.50) WITHIN GROUP (ORDER BY ne_equiv_mean_0_24h) AS median_mean,
    PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY ne_equiv_mean_0_24h) AS p75_mean,
    MAX(ne_equiv_mean_0_24h) AS max_mean,

    MIN(ne_equiv_max_0_24h) AS min_max,
    PERCENTILE_CONT(0.50) WITHIN GROUP (ORDER BY ne_equiv_max_0_24h) AS median_max,
    MAX(ne_equiv_max_0_24h) AS max_max,

    MIN(ne_equiv_auc_0_24h) AS min_auc,
    PERCENTILE_CONT(0.50) WITHIN GROUP (ORDER BY ne_equiv_auc_0_24h) AS median_auc,
    MAX(ne_equiv_auc_0_24h) AS max_auc
FROM bdmcc.vaso_ne_equiv_0_24h
WHERE ne_equiv_auc_0_24h > 0;
