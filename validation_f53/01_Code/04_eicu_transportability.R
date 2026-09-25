invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(data.table);library(sandwich);library(survival)})
logfile<-file("04_QC/eicu_transportability.log","wt");sink(logfile,split=TRUE)
e <- as.data.table(readRDS("restricted_cache/eicu_sepsis_source.rds"))
h <- as.data.table(readRDS("restricted_cache/eicu_hours.rds"))
flow<-list();record<-function(name,d)flow[[name]]<<-data.table(stage=name,n=nrow(d))
record("Local derived Sepsis-3 source",e)
e<-e[age>=18];record("Adults",e)
e<-e[icu_hadm_order==1];record("First ICU per hospitalization",e)
e<-e[unitdischargeoffset>1440 & hospitaldischargeoffset>1440];record("Still in ICU and hospital after 24h",e)
e<-e[hospitaldischargestatus %in% c("Alive","Expired")];record("Known hospital discharge vital status",e)
e<-e[icd_malignancy==0];record("No recorded malignancy",e)
e<-e[crtr_sepsis3_hr<=24];record("Local sepsis onset by landmark",e)
exp<-h[,.(mean_map=mean(map_mean),n_hour_obs=.N,n_readings=sum(readings)),by=stay_id]
e<-merge(e,exp,by="stay_id",all=FALSE);record("Any valid arterial MAP",e)
e<-e[n_hour_obs>=18];record("At least 18 observed hours",e)
e[,hospital_death:=as.integer(hospitaldischargestatus=="Expired")]
e[,strict_infection_proxy:=!is.na(antibiotic_culture_sofa_hr) & antibiotic_culture_sofa_hr<=24]
e[,gender:=factor(gender)]
e[,map_q:=relevel(cut(mean_map,c(-Inf,71.40,75.72,80.99,Inf),right=FALSE,labels=c("Q1","Q2","Q3","Q4")),"Q2")]
e[,careunit_group:=fcase(unittype %in% c("Cardiac ICU","CCU-CTICU","CSICU","CTICU"),"Cardiac ICU",
                        unittype %in% c("MICU","Med-Surg ICU"),"Medical ICU",
                        unittype %in% c("SICU","CTICU/SICU","Burn-Trauma ICU"),"Surgical/Trauma ICU",
                        unittype=="Neuro ICU","Neuro ICU/Stepdown",default="Other")]
fwrite(rbindlist(flow),"02_Results/E01_eicu_flow.csv")
fwrite(e[,.(n=.N,hospitals=uniqueN(hospitalid),patients=uniqueN(uniquepid),events=sum(hospital_death)),by=map_q],"02_Results/E02_eicu_category_events.csv")
fwrite(e[,.(n=.N),by=.(unittype,careunit_group)],"04_QC/eicu_icu_mapping.csv")
fwrite(e[,.(n=.N,events=sum(hospital_death)),by=.(strict_infection_proxy)],"04_QC/eicu_infection_proxy_audit.csv")
shared<-c("age","gender","weight","charlson","sofa","itvtn_24h_vent_tag","itvtn_24h_rrt_tag","drug_24h_sedative_tag")
fwrite(data.table(variable=shared,missing=vapply(shared,function(v)sum(is.na(e[[v]])),integer(1))),"04_QC/eicu_missing_covariates.csv")
saveRDS(e,"restricted_cache/eicu_validation_cohort.rds")

fit_logit <- function(d,covs,label,external=TRUE) {
 covs<-unique(covs)
 d<-d[complete.cases(d[,c("hospital_death","map_q",covs),with=FALSE])]
 covs<-covs[vapply(covs,function(v)uniqueN(d[[v]])>1,logical(1))]
 fml<-as.formula(paste("hospital_death ~ map_q",if(length(covs))paste("+",paste(covs,collapse=" + "))else ""))
 f<-glm(fml,data=as.data.frame(d),family=binomial())
 if(!f$converged || anyNA(coef(f)))stop("Unstable logistic model: ",label)
 v<-if(external) sandwich::vcovCL(f,cluster=list(d$hospitalid,d$uniquepid),type="HC1",fix=TRUE) else vcov(f)
 b<-coef(f);se<-sqrt(diag(v))
 # Use hospital-count degrees of freedom for the external robust interval.
 df<-if(external)uniqueN(d$hospitalid)-1 else Inf
 crit<-if(is.finite(df))qt(.975,df)else qnorm(.975)
 data.table(model=label,term=names(b),logOR=b,SE=se,OR=exp(b),lower=exp(b-crit*se),upper=exp(b+crit*se),
            p.value=2*pt(-abs(b/se),df=df),n=nrow(d),events=sum(d$hospital_death),hospitals=if(external)uniqueN(d$hospitalid)else 1L,
            endpoint="In-hospital mortality",status="exploratory transportability; derived cohort provenance incomplete")
}
results<-list()
results[["eicu_primary"]]<-fit_logit(e,shared,"eICU shared primary")
results[["eicu_context"]]<-fit_logit(e,c(shared,"careunit_group","icd_hypertension","icd_stroke","drug_24h_vaso_tag"),"eICU plus clinical context")
results[["eicu_complete24"]]<-fit_logit(e[n_hour_obs==24],shared,"eICU complete24")
results[["eicu_infection"]]<-fit_logit(e[strict_infection_proxy==TRUE],shared,"eICU antibiotic culture SOFA proxy")
results[["eicu_medical"]]<-fit_logit(e[careunit_group=="Medical ICU"],shared,"eICU medical ICU")
m<-as.data.table(readRDS("restricted_cache/mimic_enriched_preimputation.rds"))
m[,map_q:=relevel(cut(mean_observed,c(-Inf,71.40,75.72,80.99,Inf),right=FALSE,labels=c("Q1","Q2","Q3","Q4")),"Q2")]
m[,hospital_death:=hospital_expire_flag]
m[,gender:=factor(gender)]
results[["mimic_primary"]]<-fit_logit(m[sepsis_by24==TRUE & icu_hours>24],shared,"MIMIC-IV aligned early-sepsis hospital endpoint",FALSE)
results[["mimic_all"]]<-fit_logit(m,shared,"MIMIC-IV F52 cohort hospital endpoint",FALSE)
results[["mimic_complete24"]]<-fit_logit(m[sepsis_by24==TRUE & icu_hours>24 & n_hour_obs==24],shared,"MIMIC-IV aligned complete24",FALSE)
r<-rbindlist(results,fill=TRUE);fwrite(r,"02_Results/E03_harmonized_hospital_models.csv")

# Regional hospital IDs are not exported. Aggregate leave-largest-hospitals-out check.
top<-e[,.N,by=hospitalid][order(-N)][1:min(3,.N),hospitalid]
extra<-fit_logit(e[!hospitalid %in% top],shared,"eICU excluding three largest contributing hospitals")
fwrite(extra,"02_Results/E04_hospital_sensitivity.csv")
fwrite(data.table(metric=c("eICU_n","eICU_deaths","eICU_hospitals","eICU_unique_patients","repeat_patient_rows","eICU_complete24_n","eICU_proxy_n"),
                 value=c(nrow(e),sum(e$hospital_death),uniqueN(e$hospitalid),uniqueN(e$uniquepid),nrow(e)-uniqueN(e$uniquepid),sum(e$n_hour_obs==24),sum(e$strict_infection_proxy))),"04_QC/eicu_summary.csv")
print(r[term=="map_qQ4",.(model,OR,lower,upper,n,events,hospitals)])
sink();close(logfile)
