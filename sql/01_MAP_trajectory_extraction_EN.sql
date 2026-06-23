/* =====================================================================================
   Purpose:
   - Build the 24-hour landmark sepsis cohort for invasive arterial pressure analysis
   - Extract hourly invasive ABP-MAP values during the first 24 ICU hours
   - Output the design cohort, hourly MAP table, and landmark survival table

   Data source:
   - MIMIC-IV v3.1
   - Cohort source: bdmcc.bdmcc_population

   Main settings:
   - Landmark time: intime + 24 hours
   - MAP source: invasive ABP mean pressure only (itemid = 220052)
   - Hourly exposure: map_mean (primary) and map_first (sensitivity)
   - Coverage threshold: default >= 18 observed hours out of 24
   - Hourly missingness is retained for downstream multiple imputation in R

   Output tables:
   - bdmcc.cohort_design_abp
   - bdmcc.map_hourly_final_abp
   - bdmcc.base_landmark_surv_abp
===================================================================================== */

------------------------------------------------------------
-- Step 0: Drop existing output and temporary tables
------------------------------------------------------------
DROP TABLE IF EXISTS bdmcc.cohort_design_abp;
DROP TABLE IF EXISTS bdmcc.map_hourly_final_abp;
DROP TABLE IF EXISTS bdmcc.base_landmark_surv_abp;

DROP TABLE IF EXISTS tmp_param;
DROP TABLE IF EXISTS tmp_sepsis_cohort;
DROP TABLE IF EXISTS tmp_hour_grid;
DROP TABLE IF EXISTS tmp_map_raw_abp;
DROP TABLE IF EXISTS tmp_map_hourly_abp;
DROP TABLE IF EXISTS tmp_map_quality_abp;
DROP TABLE IF EXISTS tmp_keep_ids_abp;

------------------------------------------------------------
-- Step 1: Set the hourly coverage threshold
------------------------------------------------------------
CREATE TEMP TABLE tmp_param AS
SELECT
  18::int AS min_hours_required;

------------------------------------------------------------
-- Step 2: Build the design cohort and flag pre-landmark death records
------------------------------------------------------------
CREATE TEMP TABLE tmp_sepsis_cohort AS
SELECT
  p.stay_id,
  p.hadm_id,
  p.subject_id,
  p.intime,
  p.outtime,
  (p.outtime - p.intime) AS icu_los,
  (p.intime + INTERVAL '24 hours') AS landmark_time,
  p.age,
  p.crtr_sepsis3_time,
  a.deathtime AS death_time,

  CASE
    WHEN a.deathtime IS NOT NULL
     AND a.deathtime < (p.intime + INTERVAL '24 hours')
    THEN 1 ELSE 0
  END AS flag_death_lt24h,

  CASE
    WHEN a.deathtime IS NULL
      OR a.deathtime >= (p.intime + INTERVAL '24 hours')
    THEN 1 ELSE 0
  END AS landmark_eligible

FROM bdmcc.bdmcc_population p
LEFT JOIN mimiciv_hosp.admissions a
  ON a.hadm_id = p.hadm_id
WHERE
  p.crtr_sepsis3 = 1
  AND p.icu_subject_order = 1
  AND (p.outtime - p.intime) >= INTERVAL '24 hours'
  AND p.age >= 18
  AND p.icd_malignancy = 0
  AND p.icd_pregnancy = 0
;

------------------------------------------------------------
-- Step 2b: Save the full design cohort with eligibility flags
------------------------------------------------------------
CREATE TABLE bdmcc.cohort_design_abp AS
SELECT *
FROM tmp_sepsis_cohort;

------------------------------------------------------------
-- Step 3: Generate the 24 hourly bins for eligible stays
------------------------------------------------------------
CREATE TEMP TABLE tmp_hour_grid AS
SELECT
  c.stay_id,
  c.hadm_id,
  c.subject_id,
  c.intime,
  c.outtime,
  c.landmark_time,
  c.death_time,
  c.age,
  c.crtr_sepsis3_time,
  gs.hour AS hour,
  (c.intime + (gs.hour || ' hour')::interval) AS h_start,
  (c.intime + ((gs.hour + 1) || ' hour')::interval) AS h_end
FROM tmp_sepsis_cohort c
CROSS JOIN LATERAL (
  SELECT generate_series(0, 23) AS hour
) gs
WHERE c.landmark_eligible = 1;

------------------------------------------------------------
-- Step 4: Extract invasive ABP-MAP values during the first 24 ICU hours
------------------------------------------------------------
CREATE TEMP TABLE tmp_map_raw_abp AS
SELECT
  c.stay_id,
  ce.charttime,
  ce.valuenum AS map_value
FROM tmp_sepsis_cohort c
JOIN mimiciv_icu.chartevents ce
  ON ce.stay_id = c.stay_id
WHERE
  c.landmark_eligible = 1
  AND ce.charttime IS NOT NULL
  AND ce.charttime >= c.intime
  AND ce.charttime <  c.intime + INTERVAL '24 hours'
  AND ce.itemid = 220052
  AND ce.valuenum IS NOT NULL
  AND ce.valuenum BETWEEN 20 AND 200
;

------------------------------------------------------------
-- Step 5: Aggregate MAP values within each ICU hour
------------------------------------------------------------
CREATE TEMP TABLE tmp_map_hourly_abp AS
WITH joined AS (
  SELECT
    g.stay_id,
    g.hour,
    r.charttime,
    r.map_value
  FROM tmp_hour_grid g
  LEFT JOIN tmp_map_raw_abp r
    ON r.stay_id = g.stay_id
   AND r.charttime >= g.h_start
   AND r.charttime <  g.h_end
),
first_in_hour AS (
  SELECT
    stay_id,
    hour,
    map_value,
    charttime,
    ROW_NUMBER() OVER (PARTITION BY stay_id, hour ORDER BY charttime ASC) AS rn
  FROM joined
  WHERE map_value IS NOT NULL
)
SELECT
  g.stay_id,
  g.hour,
  AVG(j.map_value) AS map_mean,
  MAX(CASE WHEN f.rn = 1 THEN f.map_value END) AS map_first
FROM tmp_hour_grid g
LEFT JOIN joined j
  ON j.stay_id = g.stay_id AND j.hour = g.hour
LEFT JOIN first_in_hour f
  ON f.stay_id = g.stay_id AND f.hour = g.hour
GROUP BY
  g.stay_id, g.hour
ORDER BY
  g.stay_id, g.hour
;

------------------------------------------------------------
-- Step 6: Count observed hourly MAP values per stay
------------------------------------------------------------
CREATE TEMP TABLE tmp_map_quality_abp AS
SELECT
  stay_id,
  COUNT(*) FILTER (WHERE map_mean IS NOT NULL) AS n_hours_mean
FROM tmp_map_hourly_abp
GROUP BY stay_id;

------------------------------------------------------------
-- Step 7: Keep stays meeting the hourly coverage threshold
------------------------------------------------------------
CREATE TEMP TABLE tmp_keep_ids_abp AS
SELECT
  q.stay_id
FROM tmp_map_quality_abp q
CROSS JOIN tmp_param p
WHERE
  q.n_hours_mean >= p.min_hours_required;

------------------------------------------------------------
-- Step 8: Save the final hourly MAP table
------------------------------------------------------------
CREATE TABLE bdmcc.map_hourly_final_abp AS
SELECT
  m.stay_id,
  m.hour,
  m.map_mean,
  m.map_first
FROM tmp_map_hourly_abp m
JOIN tmp_keep_ids_abp k
  ON k.stay_id = m.stay_id
ORDER BY
  m.stay_id, m.hour;

------------------------------------------------------------
-- Step 9: Save the landmark survival table
------------------------------------------------------------
CREATE TABLE bdmcc.base_landmark_surv_abp AS
SELECT
  c.stay_id,
  c.subject_id,
  c.hadm_id,
  c.intime,
  c.outtime,
  c.landmark_time,
  c.death_time,
  c.age,
  c.crtr_sepsis3_time,
  c.flag_death_lt24h,
  c.landmark_eligible,
  q.n_hours_mean,

  CASE
    WHEN c.death_time IS NOT NULL
     AND c.death_time >= c.landmark_time
     AND c.death_time <= c.outtime
    THEN 1 ELSE 0
  END AS event_icu_death_after24h,

  EXTRACT(EPOCH FROM (
    (
      CASE
        WHEN c.death_time IS NOT NULL
         AND c.death_time >= c.landmark_time
         AND c.death_time <= c.outtime
        THEN c.death_time
        ELSE c.outtime
      END
    ) - c.landmark_time
  )) / 86400.0 AS time_from_24h_days

FROM tmp_sepsis_cohort c
JOIN tmp_keep_ids_abp k
  ON k.stay_id = c.stay_id
LEFT JOIN tmp_map_quality_abp q
  ON q.stay_id = c.stay_id
WHERE
  c.landmark_eligible = 1
;

------------------------------------------------------------
-- Step 10: Post-run quality checks
------------------------------------------------------------

-- 1) Confirm that the output tables were created
SELECT to_regclass('bdmcc.cohort_design_abp')       AS design_table,
       to_regclass('bdmcc.map_hourly_final_abp')    AS map_table,
       to_regclass('bdmcc.base_landmark_surv_abp')  AS surv_table;

-- 2) Summarize the design cohort and landmark eligibility flags
SELECT
  COUNT(*) AS n_design,
  SUM(flag_death_lt24h) AS n_flag_death_lt24h,
  SUM(landmark_eligible) AS n_landmark_eligible
FROM bdmcc.cohort_design_abp;

-- 3) Count unique stays in the final analysis cohort
SELECT COUNT(DISTINCT stay_id) AS n_analysis
FROM bdmcc.base_landmark_surv_abp;

-- 4) Confirm that no ineligible stays remain in the landmark table
SELECT COUNT(*) AS n_bad_landmark_eligible
FROM bdmcc.base_landmark_surv_abp
WHERE landmark_eligible = 0;

-- 5) Check for negative follow-up time
SELECT COUNT(*) AS n_negative_time
FROM bdmcc.base_landmark_surv_abp
WHERE time_from_24h_days < 0;

-- 6) Check the number of hourly rows per stay
SELECT stay_id, COUNT(*) AS n_rows
FROM bdmcc.map_hourly_final_abp
GROUP BY stay_id
ORDER BY n_rows DESC
LIMIT 10;

-- 7) Summarize hourly coverage under common thresholds
SELECT
  COUNT(*) AS n_total,
  SUM(CASE WHEN n_hours_mean >= 16 THEN 1 ELSE 0 END) AS n_ge_16,
  SUM(CASE WHEN n_hours_mean >= 18 THEN 1 ELSE 0 END) AS n_ge_18,
  SUM(CASE WHEN n_hours_mean >= 16 AND n_hours_mean < 18 THEN 1 ELSE 0 END) AS n_16_17
FROM bdmcc.base_landmark_surv_abp;
