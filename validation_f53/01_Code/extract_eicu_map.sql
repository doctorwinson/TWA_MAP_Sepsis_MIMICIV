WITH eligible AS MATERIALIZED (
  SELECT p.stay_id
  FROM bdmcc.bdmcc_population p
  JOIN eicu_icu.patient u ON u.patientunitstayid=p.stay_id
  WHERE p.crtr_sepsis3=1 AND p.age>=18 AND p.icu_hadm_order=1
    AND u.unitdischargeoffset>1440 AND u.hospitaldischargeoffset>1440
    AND u.hospitaldischargestatus IN ('Alive','Expired')
    AND p.icd_malignancy=0 AND p.crtr_sepsis3_hr<=24
)
SELECT v.patientunitstayid AS stay_id,
       floor(v.observationoffset/60.0)::integer AS hour,
       avg(v.systemicmean)::double precision AS map_mean,
       count(*) AS readings
FROM eligible e JOIN eicu_icu.vitalperiodic v ON v.patientunitstayid=e.stay_id
WHERE v.observationoffset>=0 AND v.observationoffset<1440
  AND v.systemicmean BETWEEN 20 AND 200
GROUP BY v.patientunitstayid,floor(v.observationoffset/60.0)::integer;
