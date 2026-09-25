invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(DBI);library(odbc);library(data.table)})
setDTthreads(2)
for(path in c("restricted_cache","02_Results","04_QC"))dir.create(path,showWarnings=FALSE,recursive=TRUE)
con<-dbConnect(odbc(),dsn=Sys.getenv("MIMIC_DSN","mimic4_v31"),timeout=15)
dbExecute(con,"SET default_transaction_read_only=on")
dbExecute(con,"SET statement_timeout='20min'")
q<-function(s)as.data.table(dbGetQuery(con,s))
d<-q(paste(readLines("01_Code/19_rebuild_source.sql"),collapse="\n"))
stopifnot(!anyDuplicated(d$stay_id),!anyNA(d$crtr_sepsis3_time))
flow<-data.table(stage="Sepsis-3 stays",n=nrow(d))
restrict<-function(x,label){flow<<-rbind(flow,data.table(stage=label,n=nrow(x)));x}
d<-restrict(d[first_rank==1],"First recorded ICU stay")
d[,icu_hours:=as.numeric(difftime(outtime,intime,units="hours"))]
d<-restrict(d[age>=18 & icu_hours>24],"Adult and still ICU after hour24")
d<-restrict(d[icd_malignancy==0 & icd_pregnancy==0],"No coded malignancy or pregnancy")
d[,onset_hours:=as.numeric(difftime(crtr_sepsis3_time,intime,units="hours"))]
d<-restrict(d[onset_hours<=24],"Sepsis onset by hour24")
d[,death_time:=deathtime]
d[is.na(death_time),death_time:=as.POSIXct(dod,tz="UTC")]
d[,death_hours:=as.numeric(difftime(death_time,intime,units="hours"))]
d<-restrict(d[is.na(death_hours)|death_hours>24],"Alive at hour24")
d[,event_lm:=as.integer(!is.na(death_hours)&death_hours<=720)]
d[,time_lm_days:=fifelse(event_lm==1,(death_hours-24)/24,29)]
d[,both_pair_by24:=!is.na(antibiotic_time)&!is.na(culture_time)&antibiotic_time<=intime+86400&culture_time<=intime+86400]
ids<-paste(d$stay_id,collapse=",")
cat("Extracting invasive MAP for",nrow(d),"stays\n")
h<-q(sprintf("SELECT i.stay_id,floor(extract(epoch FROM(ce.charttime-i.intime))/3600)::int AS hr,avg(ce.valuenum) AS value FROM mimiciv_icu.icustays i JOIN mimiciv_icu.chartevents ce USING(stay_id) WHERE i.stay_id IN (%s) AND ce.itemid=220052 AND ce.charttime>=i.intime AND ce.charttime<i.intime+interval '24 hours' AND ce.valuenum BETWEEN 20 AND 200 GROUP BY i.stay_id,hr",ids))
setnames(h,"hr","hour")
saveRDS(h,"restricted_cache/F55_raw_hourly_MAP.rds")
hours<-sprintf("h%02d",0:23)
h[,hour:=sprintf("h%02d",hour)]
wide<-dcast(h,stay_id~hour,value.var="value")
d<-merge(d,wide,by="stay_id",all.x=TRUE)
d[,n_hour_obs:=rowSums(!is.na(.SD)),.SDcols=hours]
d<-restrict(d[n_hour_obs>=18],"At least18 observed hours")
d[,mean_observed:=rowMeans(.SD,na.rm=TRUE),.SDcols=hours]
ids<-paste(d$stay_id,collapse=",")
cat("Extracting interventions for",nrow(d),"stays\n")
v<-q(sprintf("SELECT i.stay_id,max(CASE WHEN v.ventilation_status='InvasiveVent' THEN 1 ELSE 0 END) invasive,max(CASE WHEN v.ventilation_status IN ('InvasiveVent','Tracheostomy') THEN 1 ELSE 0 END) invasive_or_trach FROM mimiciv_icu.icustays i LEFT JOIN mimiciv_derived.ventilation v ON v.stay_id=i.stay_id AND v.starttime<i.intime+interval '24 hours' AND v.endtime>i.intime WHERE i.stay_id IN (%s) GROUP BY i.stay_id",ids))
d<-merge(d,v,by="stay_id",all.x=TRUE)
d[,itvtn_24h_vent_tag:=invasive]
# Tracheostomy alone does not establish mechanical ventilation.
sed<-c(221385,221623,221668,221712,221744,221833,222168,225150,225154,225156,225942,225972,229420)
ie<-q(sprintf("SELECT ie.stay_id,ie.itemid,ie.starttime,ie.endtime,ie.amount,ie.amountuom,ie.rate,ie.rateuom,ie.patientweight,ie.statusdescription FROM mimiciv_icu.inputevents ie JOIN mimiciv_icu.icustays i USING(stay_id) WHERE ie.stay_id IN (%s) AND ie.itemid IN (%s,221906,221289,221662,221749,222315) AND ie.starttime<i.intime+interval '24 hours' AND coalesce(ie.endtime,ie.starttime)>=i.intime",ids,paste(sed,collapse=",")))
saveRDS(ie,"restricted_cache/F55_drug_events.rds")
ie<-merge(ie,d[,.(stay_id,intime)],by="stay_id")
ie<-ie[is.na(statusdescription)|statusdescription!="Rewritten"]
ie[,duration_h:=pmax(0,as.numeric(difftime(pmin(endtime,intime+86400),pmax(starttime,intime),units="hours")))]
valid_drug<-ie[itemid %in% sed & (amount>0|rate>0) & (duration_h>0 | (starttime>=intime & starttime<intime+86400)),unique(stay_id)]
d[,drug_24h_sedative_tag:=as.integer(stay_id %in% valid_drug)]
vp<-ie[itemid %in% c(221906,221289,221662,221749,222315)&rate>0&duration_h>0]
vp[,unit:=tolower(rateuom)]
vp[,std:=fcase(itemid==222315 & unit %in% c("units/min","unit/min"),rate,
 itemid==222315 & unit %in% c("units/hour","unit/hour"),rate/60,
 itemid!=222315 & unit %in% c("mcg/kg/min","ug/kg/min"),rate,
 itemid!=222315 & unit %in% c("mcg/min","ug/min") & patientweight>0,rate/patientweight,
 default=NA_real_)]
fwrite(vp[,.(records=.N,unhandled=sum(is.na(std))),by=.(itemid,unit)],"04_QC/F55_vasopressor_unit_audit.csv")
stopifnot(!anyNA(vp$std))
vp[,factor:=fcase(itemid==221662,.01,itemid==221749,.1,itemid==222315,2.5,default=1)]
vs<-vp[,.(ne_equiv_mean_0_24h=sum(std*factor*duration_h)/24),by=stay_id]
d<-merge(d,vs,by="stay_id",all.x=TRUE)
d[is.na(ne_equiv_mean_0_24h),ne_equiv_mean_0_24h:=0]
d[,ne_log1p:=log1p(ne_equiv_mean_0_24h)]
mon<-q(sprintf("SELECT i.stay_id,max(CASE WHEN ce.itemid=224322 AND ce.valuenum>0 THEN 1 ELSE 0 END) iabp,max(CASE WHEN ce.itemid BETWEEN 228176 AND 228185 AND ce.valuenum IS NOT NULL THEN 1 ELSE 0 END) picco,max(CASE WHEN ce.itemid IN (228368,228369,228371,228374,228375,228376,228378,228379,228380,228381,228382) AND ce.valuenum IS NOT NULL THEN 1 ELSE 0 END) nicom FROM mimiciv_icu.icustays i LEFT JOIN mimiciv_icu.chartevents ce ON ce.stay_id=i.stay_id AND ce.charttime>=i.intime AND ce.charttime<i.intime+interval '24 hours' AND (ce.itemid=224322 OR ce.itemid BETWEEN 228176 AND 228185 OR ce.itemid IN (228368,228369,228371,228374,228375,228376,228378,228379,228380,228381,228382)) WHERE i.stay_id IN (%s) GROUP BY i.stay_id",ids))
pr<-q(sprintf("SELECT DISTINCT i.stay_id FROM mimiciv_icu.icustays i JOIN mimiciv_icu.procedureevents pe USING(stay_id) WHERE i.stay_id IN (%s) AND pe.itemid=224272 AND pe.starttime<i.intime+interval '24 hours' AND pe.endtime>i.intime",ids))
d<-merge(d,mon,by="stay_id",all.x=TRUE)
d[,itvtn_24h_iabp_tag:=as.integer(iabp==1|stay_id %in% pr$stay_id)]
d[,advanced_monitor:=as.integer(picco==1|nicom==1)]
d[,cardio_proc_device:=as.integer(itvtn_24h_iabp_tag==1|itvtn_24h_cabg_tag==1|itvtn_24h_pci_tag==1)]
lac<-q(sprintf("SELECT DISTINCT ON(i.stay_id) i.stay_id,l.valuenum AS lab_fst_24h_lactate_first FROM mimiciv_icu.icustays i JOIN mimiciv_hosp.labevents l ON l.hadm_id=i.hadm_id AND l.charttime>=i.intime AND l.charttime<i.intime+interval '24 hours' WHERE i.stay_id IN (%s) AND l.itemid=50813 AND l.valuenum>0 ORDER BY i.stay_id,l.charttime,l.labevent_id",ids))
d<-merge(d,lac,by="stay_id",all.x=TRUE)
d[,careunit_group:=factor(fcase(first_careunit %in% c("Cardiac Vascular Intensive Care Unit (CVICU)","Coronary Care Unit (CCU)"),"Cardiac ICU",first_careunit %in% c("Medical Intensive Care Unit (MICU)","Medical/Surgical Intensive Care Unit (MICU/SICU)"),"Medical ICU",first_careunit %in% c("Surgical Intensive Care Unit (SICU)","Trauma SICU (TSICU)","PACU","Surgery/Vascular/Intermediate"),"Surgical/Trauma ICU",first_careunit %in% c("Neuro Intermediate","Neuro Surgical Intensive Care Unit (Neuro SICU)","Neuro Stepdown","Neurology"),"Neuro ICU/Stepdown",default="Other"),levels=c("Medical ICU","Cardiac ICU","Surgical/Trauma ICU","Neuro ICU/Stepdown","Other"))]
d[,admission_group:=factor(fcase(admission_type %in% c("EW EMER.","DIRECT EMER."),"Emergency/Direct emergent",admission_type=="URGENT","Urgent",admission_type %in% c("ELECTIVE","SURGICAL SAME DAY ADMISSION"),"Elective/Same-day surgery",admission_type %in% c("EU OBSERVATION","OBSERVATION ADMIT","DIRECT OBSERVATION","AMBULATORY OBSERVATION"),"Observation",default="Other"),levels=c("Emergency/Direct emergent","Urgent","Elective/Same-day surgery","Observation","Other"))]
d[,service_group:=factor(fcase(service_at_icu_entry %in% c("NMED","NSURG"),"Neuro service",service_at_icu_entry %in% c("CSURG","SURG","TSURG","VSURG","ORTHO","PSURG","ENT","GU","GYN"),"Surgical service",service_at_icu_entry %in% c("CMED","MED","OMED"),"Medical service",default="Other/Unavailable"),levels=c("Medical service","Surgical service","Neuro service","Other/Unavailable"))]
d[,year_offset:=as.integer(format(intime,"%Y"))-anchor_year]
d[,year_lower:=as.integer(substr(anchor_year_group,1,4))+year_offset]
d[,year_upper:=as.integer(substr(anchor_year_group,8,11))+year_offset]
d[,period:=factor(fcase(year_upper<=2016,"Earlier <=2016",year_lower>=2017,"Later >=2017",default="Boundary interval"),levels=c("Earlier <=2016","Later >=2017","Boundary interval"))]
primary<-c("age","gender","weight","charlson","sofa","itvtn_24h_vent_tag","itvtn_24h_rrt_tag","drug_24h_sedative_tag")
fwrite(d[,lapply(.SD,function(x)sum(is.na(x))),.SDcols=primary],"04_QC/primary_missingness_before_exclusion.csv")
saveRDS(d,"restricted_cache/source_enriched.rds")
saveRDS(list(primary=primary,hours=hours),"restricted_cache/pre_spec.rds")
saveRDS(d[,.(stay_id,both_pair_by24)],"restricted_cache/sepsis_source_verification.rds")
final<-restrict(d[complete.cases(d[,..primary])],"Complete primary covariates")
stopifnot(!anyDuplicated(final$subject_id),all(final$time_lm_days>0 & final$time_lm_days<=29))
fwrite(flow,"02_Results/F55_flow.csv")
fwrite(data.table(stage=c("Reconstructed observed18 source","Complete primary covariates"),n=c(nrow(d),nrow(final))),"02_Results/cohort_counts.csv")
legacy_path<-Sys.getenv("F55_LEGACY_AUDIT_RDS",file.path(dirname(getwd()),"增强转投_F54_20260925_landmark_corrected/restricted_cache/eligible_primary.rds"))
if(file.exists(legacy_path)) {
legacy<-as.data.table(readRDS(legacy_path))
comp<-merge(legacy,final,by="stay_id",suffixes=c("_F54","_F55"))
vars<-c("age","weight","icd_stroke","icd_hypertension","itvtn_24h_vent_tag","drug_24h_sedative_tag","itvtn_24h_cabg_tag","itvtn_24h_pci_tag","itvtn_24h_iabp_tag","advanced_monitor","ne_equiv_mean_0_24h")
audit<-rbindlist(lapply(vars,function(v){a<-comp[[paste0(v,"_F54")]];b<-comp[[paste0(v,"_F55")]];data.table(variable=v,paired=sum(!is.na(a)&!is.na(b)),different=sum(abs(a-b)>1e-6,na.rm=TRUE))}))
fwrite(audit,"04_QC/F55_vs_F54_variable_changes.csv")
fwrite(data.table(metric=c("F54_n","F55_n","shared","F54_only","F55_only"),n=c(nrow(legacy),nrow(final),nrow(comp),sum(!legacy$stay_id %in% final$stay_id),sum(!final$stay_id %in% legacy$stay_id))),"04_QC/F55_vs_F54_cohort_changes.csv")
print(audit)
}
print(flow)
dbDisconnect(con)
