WITH candidates AS (
 SELECT c.stay_id,c.intime
 FROM bdmcc.cohort_design_abp c
 JOIN bdmcc.bdmcc_population p USING(stay_id)
 WHERE c.landmark_eligible=1
   AND p.crtr_sepsis3_time<=c.intime+INTERVAL '24 hours'
)
SELECT c.stay_id,
 floor(extract(epoch FROM (ce.charttime-c.intime))/3600)::integer AS hour,
 avg(ce.valuenum) AS independently_extracted_map
FROM candidates c
JOIN mimiciv_icu.chartevents ce USING(stay_id)
WHERE ce.itemid=220052
 AND ce.charttime>=c.intime
 AND ce.charttime<c.intime+INTERVAL '24 hours'
 AND ce.valuenum BETWEEN 20 AND 200
GROUP BY c.stay_id,hour;
