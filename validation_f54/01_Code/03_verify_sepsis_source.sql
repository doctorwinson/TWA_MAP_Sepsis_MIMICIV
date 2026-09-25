SELECT p.stay_id,
       p.crtr_sepsis3_time AS source_onset,
       p.crtr_sepsis3_suspected_infection_time AS source_infection_time,
       p.crtr_sepsis3_sofa_time AS source_sofa_time,
       s.suspected_infection_time, s.sofa_time, s.sofa_score,
       s.antibiotic_time, s.culture_time, s.sepsis3,
       p.intime, p.day30_outcome, p.day30_los,
       a.deathtime, h.dod
FROM bdmcc.bdmcc_population p
LEFT JOIN mimiciv_derived.sepsis3 s USING(stay_id)
LEFT JOIN mimiciv_hosp.admissions a USING(hadm_id)
LEFT JOIN mimiciv_hosp.patients h ON h.subject_id=p.subject_id
WHERE p.crtr_sepsis3=1;
