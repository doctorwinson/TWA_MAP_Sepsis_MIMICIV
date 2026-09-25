-- PostgreSQL; read-only. No bdmcc dependency.
-- Hospital diagnosis/procedure codes describe admission-level context, not
-- diagnoses or procedures known to precede the 24-hour landmark.
WITH ranked AS (
 SELECT i.*,row_number() OVER(PARTITION BY subject_id ORDER BY intime,stay_id) first_rank
 FROM mimiciv_icu.icustays i
), dx AS (
 SELECT hadm_id,
 max(CASE WHEN (icd_version=9 AND (substring(icd_code,1,3) BETWEEN '140' AND '208' OR substring(icd_code,1,4) IN ('2090','2091','2092','2093','2097'))) OR (icd_version=10 AND icd_code LIKE 'C%') THEN 1 ELSE 0 END) icd_malignancy,
 max(CASE WHEN (icd_version=9 AND (substring(icd_code,1,3) BETWEEN '630' AND '679' OR substring(icd_code,1,3) IN ('V22','V23','V27'))) OR (icd_version=10 AND (icd_code LIKE 'O%' OR substring(icd_code,1,3) IN ('Z33','Z34','Z35','Z37'))) THEN 1 ELSE 0 END) icd_pregnancy,
 max(CASE WHEN (icd_version=9 AND substring(icd_code,1,3) BETWEEN '430' AND '438') OR (icd_version=10 AND substring(icd_code,1,3) BETWEEN 'I60' AND 'I69') THEN 1 ELSE 0 END) icd_stroke,
 max(CASE WHEN (icd_version=9 AND substring(icd_code,1,3) BETWEEN '401' AND '405') OR (icd_version=10 AND substring(icd_code,1,3) BETWEEN 'I10' AND 'I16') THEN 1 ELSE 0 END) icd_hypertension,
 max(CASE WHEN (icd_version=9 AND substring(icd_code,1,3) BETWEEN '410' AND '414') OR (icd_version=10 AND substring(icd_code,1,3) BETWEEN 'I20' AND 'I25') THEN 1 ELSE 0 END) icd_cad
 FROM mimiciv_hosp.diagnoses_icd GROUP BY hadm_id
), proc AS (
 SELECT hadm_id,
 max(CASE WHEN (icd_version=9 AND icd_code LIKE '361%') OR (icd_version=10 AND substring(icd_code,1,4) IN ('0210','0211','0212','0213')) THEN 1 ELSE 0 END) cabg_code,
 max(CASE WHEN (icd_version=9 AND icd_code IN ('0066','3601','3602','3605','3606','3607')) OR (icd_version=10 AND substring(icd_code,1,4) IN ('0270','0271','0272','0273')) THEN 1 ELSE 0 END) pci_code
 FROM mimiciv_hosp.procedures_icd GROUP BY hadm_id
)
SELECT i.subject_id,i.hadm_id,i.stay_id,i.intime,i.outtime,i.first_rank,i.first_careunit,
 a.admittime,a.admission_type,a.hospital_expire_flag,a.deathtime,p.dod,p.gender,
 p.anchor_age,p.anchor_year,p.anchor_year_group,
 p.anchor_age + extract(year FROM i.intime)::int-p.anchor_year AS age,
 s.sofa_time,s.suspected_infection_time,s.antibiotic_time,s.culture_time,
 greatest(s.sofa_time,s.suspected_infection_time) AS crtr_sepsis3_time,
 coalesce(dx.icd_malignancy,0) icd_malignancy,coalesce(dx.icd_pregnancy,0) icd_pregnancy,
 coalesce(dx.icd_stroke,0) icd_stroke,coalesce(dx.icd_hypertension,0) icd_hypertension,coalesce(dx.icd_cad,0) icd_cad,
 coalesce(proc.cabg_code,0) itvtn_24h_cabg_tag,coalesce(proc.pci_code,0) itvtn_24h_pci_tag,
 w.weight,w.weight_admit,c.charlson_comorbidity_index AS charlson,fs.sofa,
 fs.respiration sofa_respiration,fs.coagulation sofa_coagulation,fs.liver sofa_liver,
 fs.cardiovascular sofa_cardiovascular,fs.cns sofa_cns,fs.renal sofa_renal,
 coalesce(r.dialysis_active,0) itvtn_24h_rrt_tag,sp.sapsii,
 sv.curr_service service_at_icu_entry
FROM ranked i
JOIN mimiciv_derived.sepsis3 s USING(stay_id)
JOIN mimiciv_hosp.patients p ON p.subject_id=i.subject_id
JOIN mimiciv_hosp.admissions a ON a.hadm_id=i.hadm_id
LEFT JOIN dx ON dx.hadm_id=i.hadm_id
LEFT JOIN proc ON proc.hadm_id=i.hadm_id
LEFT JOIN mimiciv_derived.first_day_weight w USING(stay_id)
LEFT JOIN mimiciv_derived.charlson c ON c.hadm_id=i.hadm_id
LEFT JOIN mimiciv_derived.first_day_sofa fs USING(stay_id)
LEFT JOIN mimiciv_derived.first_day_rrt r USING(stay_id)
LEFT JOIN mimiciv_derived.sapsii sp USING(stay_id)
LEFT JOIN LATERAL (SELECT curr_service FROM mimiciv_hosp.services x WHERE x.hadm_id=i.hadm_id AND x.transfertime<=i.intime ORDER BY transfertime DESC,curr_service LIMIT 1) sv ON true
WHERE s.sepsis3=true;
