-- Local derived definitions are provisional until their upstream SQL is recovered.
-- Never infer complete 30-day vital status from hospital discharge status.
SELECT p.stay_id,p.subject_id,p.hadm_id,p.age,p.gender,p.weight,p.charlson,p.sofa,
       p.icu_hadm_order,p.icu_los,p.hos_los,p.icd_malignancy,
       p.crtr_sepsis3,p.crtr_sepsis3_hr,
       p.itvtn_24h_vent_tag,p.itvtn_24h_rrt_tag,p.drug_24h_sedative_tag,
       p.drug_24h_vaso_tag,p.icd_hypertension,p.icd_stroke,
       u.hospitalid,u.uniquepid,u.patienthealthsystemstayid,u.unittype,
       u.unitadmitSource,u.unitdischargeoffset,u.hospitaldischargeoffset,
       u.hospitaldischargestatus,u.unitdischargestatus,
       s.antibiotic_culture_sofa_hr,s.infection_diagnosis_sofa_hr,
       s.sepsis_diagnosis_offset_hr,s.severe_sepsis_onset_hr,
       s.suspected_infection_hr,s.organ_dysfunction_hr
FROM bdmcc.bdmcc_population p
JOIN eicu_icu.patient u ON u.patientunitstayid=p.stay_id
LEFT JOIN bdmcc.sepsis3 s ON s.patientunitstayid=p.stay_id
WHERE p.crtr_sepsis3=1;
