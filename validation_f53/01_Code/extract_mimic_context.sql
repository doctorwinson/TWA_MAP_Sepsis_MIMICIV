SELECT p.stay_id, pa.anchor_year, pa.anchor_year_group,
       i.first_careunit, a.admission_type, a.hospital_expire_flag,
       a.deathtime, a.dischtime AS hospital_dischtime,
       fs.sofa AS sofa_domains_total,
       fs.respiration AS sofa_respiration, fs.coagulation AS sofa_coagulation,
       fs.liver AS sofa_liver, fs.cardiovascular AS sofa_cardiovascular,
       fs.cns AS sofa_cns, fs.renal AS sofa_renal,
       s.curr_service AS service_at_icu_entry,
       s.transfertime AS service_timestamp
FROM bdmcc.base_landmark_surv_abp b
JOIN bdmcc.bdmcc_population p ON p.stay_id=b.stay_id
JOIN mimiciv_hosp.patients pa ON pa.subject_id=p.subject_id
JOIN mimiciv_icu.icustays i ON i.stay_id=p.stay_id
JOIN mimiciv_hosp.admissions a ON a.hadm_id=p.hadm_id
LEFT JOIN mimiciv_derived.first_day_sofa fs ON fs.stay_id=p.stay_id
LEFT JOIN LATERAL (
  SELECT curr_service, transfertime
  FROM mimiciv_hosp.services sv
  WHERE sv.hadm_id=p.hadm_id AND sv.transfertime<=p.intime
  ORDER BY transfertime DESC, curr_service ASC LIMIT 1
) s ON true;
