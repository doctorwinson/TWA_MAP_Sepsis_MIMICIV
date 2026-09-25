invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(data.table);library(mice);library(survival);library(ggplot2)})
options(width=140)
dir.create("02_Results",showWarnings=FALSE)
logfile <- file("04_QC/phenotype_validation.log","wt");sink(logfile,split=TRUE)
cat("Started",format(Sys.time()),"\n")
src <- Sys.getenv("F53_SOURCE_ROOT",dirname(normalizePath(getwd(),winslash="/")))
raw <- as.data.table(readRDS(file.path(src,"outputs_abp_map/cache/cache_main_ge18h_Model2_SOFA_minimal.rds"))$dat_lm18)
context <- readRDS("restricted_cache/mimic_context.rds")
mi <- readRDS("restricted_cache/mice40_full.rds")
spec <- readRDS("restricted_cache/model_spec.rds")
stopifnot(mi$imp$m==40L)
primary <- spec$covars
domains <- paste0("sofa_",c("respiration","coagulation","liver","cardiovascular","cns","renal"))
hour_cols <- spec$hour_cols

classify_icu <- function(x) {
  out <- rep("Other",length(x))
  out[x %in% c("Cardiac Vascular Intensive Care Unit (CVICU)","Coronary Care Unit (CCU)")] <- "Cardiac ICU"
  out[x %in% c("Medical Intensive Care Unit (MICU)","Medical/Surgical Intensive Care Unit (MICU/SICU)")] <- "Medical ICU"
  out[x %in% c("Surgical Intensive Care Unit (SICU)","Trauma SICU (TSICU)","PACU","Surgery/Vascular/Intermediate")] <- "Surgical/Trauma ICU"
  out[x %in% c("Neuro Intermediate","Neuro Surgical Intensive Care Unit (Neuro SICU)","Neuro Stepdown","Neurology")] <- "Neuro ICU/Stepdown"
  factor(out,levels=c("Medical ICU","Cardiac ICU","Surgical/Trauma ICU","Neuro ICU/Stepdown","Other"))
}
classify_admission <- function(x) {
  out <- rep("Other",length(x))
  out[x %in% c("EW EMER.","DIRECT EMER.")] <- "Emergency/Direct emergent"
  out[x=="URGENT"] <- "Urgent"
  out[x %in% c("ELECTIVE","SURGICAL SAME DAY ADMISSION")] <- "Elective/Same-day surgery"
  out[x %in% c("EU OBSERVATION","OBSERVATION ADMIT","DIRECT OBSERVATION","AMBULATORY OBSERVATION")] <- "Observation"
  factor(out,levels=c("Emergency/Direct emergent","Urgent","Elective/Same-day surgery","Observation","Other"))
}
dat <- merge(raw,context,by="stay_id",all.x=TRUE,sort=FALSE)
stopifnot(nrow(dat)==nrow(raw),!anyDuplicated(dat$subject_id))
dat[,careunit_group:=classify_icu(first_careunit)]
dat[,admission_group:=classify_admission(admission_type)]
dat[,cardio_proc_device:=as.integer(itvtn_24h_cabg_tag==1 | itvtn_24h_pci_tag==1 | itvtn_24h_iabp_tag==1)]
dat[,advanced_monitor:=as.integer(itvtn_24h_picco_tag==1 | itvtn_24h_nicom_tag==1)]
dat[,ne_log1p:=log1p(ne_equiv_mean_0_24h)]
dat[,sofa_noncv:=rowSums(.SD),.SDcols=setdiff(domains,"sofa_cardiovascular")]
dat[,sofa_noncns:=rowSums(.SD),.SDcols=setdiff(domains,"sofa_cns")]
dat[,sofa_sum:=rowSums(.SD),.SDcols=domains]
dat[,onset_hours:=as.numeric(difftime(crtr_sepsis3_time,intime,units="hours"))]
dat[,sepsis_by24:=!is.na(onset_hours) & onset_hours<=24]
dat[,icu_hours:=as.numeric(difftime(outtime,intime,units="hours"))]
dat[,year_offset:=as.integer(format(intime,"%Y"))-anchor_year]
dat[,year_lower:=as.integer(substr(anchor_year_group,1,4))+year_offset]
dat[,year_upper:=as.integer(substr(anchor_year_group,8,11))+year_offset]
dat[,period:=fcase(year_upper<=2016,"Earlier <=2016",year_lower>=2017,"Later >=2017",default="Boundary interval")]
dat[,period:=factor(period,levels=c("Earlier <=2016","Later >=2017","Boundary interval"))]
dat[,service_group:=fcase(service_at_icu_entry %in% c("NMED","NSURG"),"Neuro service",
                        service_at_icu_entry %in% c("CSURG","SURG","TSURG","VSURG","ORTHO","PSURG","ENT","GU","GYN"),"Surgical service",
                        service_at_icu_entry %in% c("CMED","MED","OMED"),"Medical service",default="Other/Unavailable")]
dat[,service_group:=factor(service_group,levels=c("Medical service","Surgical service","Neuro service","Other/Unavailable"))]
dat[,mean_observed:=rowMeans(.SD,na.rm=TRUE),.SDcols=hour_cols]

audit <- data.table(metric=c("source_n","post_MI_n","sepsis_after_24h","sepsis_timing_missing","ICU_LOS_exactly_24h","SOFA_domain_sum_differs","service_missing","source_endpoint_disagreement","followup_disagreement"),
 value=c(nrow(dat),nrow(mi$imp$data),sum(dat$onset_hours>24,na.rm=TRUE),sum(is.na(dat$onset_hours)),sum(dat$icu_hours==24,na.rm=TRUE),
 sum(dat$sofa!=dat$sofa_sum,na.rm=TRUE),sum(is.na(dat$service_at_icu_entry)),
 sum(dat$event_lm!=as.integer(dat$day30_outcome==1 & dat$day30_los>24 & dat$day30_los<=720),na.rm=TRUE),
 sum(abs(dat$time_lm_days-(pmin(dat$day30_los,720)-24)/24)>1e-5,na.rm=TRUE)))
fwrite(audit,"04_QC/cohort_audit.csv")
fwrite(dat[,.(n=.N,events=sum(event_lm)),by=.(period)],"02_Results/V01_temporal_counts.csv")
fwrite(dat[,.(n=.N),by=.(first_careunit,careunit_group)],"04_QC/careunit_mapping_counts.csv")
fwrite(dat[,.(n=.N),by=.(service_at_icu_entry,service_group)],"04_QC/service_mapping_counts.csv")
context_keep <- c("subject_id","hospital_expire_flag","sepsis_by24","onset_hours","icu_hours","period","service_group","careunit_group","admission_group",
 "cardio_proc_device","advanced_monitor","ne_log1p","ne_equiv_mean_0_24h","icd_hypertension","icd_stroke","drug_24h_vaso_tag",domains,"sofa_noncv","sofa_noncns","sofa_sum","n_hour_obs")
extra_names <- setdiff(context_keep,names(mi$imp$data))
extras <- dat[,c("stay_id",extra_names),with=FALSE]
ddlist <- lapply(1:40,function(i){
  dd <- as.data.table(complete(mi$imp,i))
  dd <- merge(dd,extras,by="stay_id",all.x=TRUE,sort=FALSE)
  dd[,mean_map:=rowMeans(.SD),.SDcols=hour_cols]
  dd[,map_q:=relevel(cut(mean_map,c(-Inf,spec$cutpoints[2:4],Inf),labels=c("Q1","Q2","Q3","Q4")),"Q2")]
  dd[,gender:=factor(gender)]
  dd
})
saveRDS(dat,"restricted_cache/mimic_enriched_preimputation.rds")
saveRDS(ddlist,"restricted_cache/mimic_m40_enriched.rds")

get_formula <- function(covs,exposure="map_q") as.formula(paste("Surv(time_lm_days,event_lm) ~",paste(c(exposure,covs),collapse=" + ")))
fit_one <- function(dd,covs,exposure="map_q") {
  cc <- dd[complete.cases(dd[,unique(c("time_lm_days","event_lm","map_q",covs)),with=FALSE])]
  covs <- covs[vapply(covs,function(x) uniqueN(cc[[x]])>1L,logical(1))]
  cc <- droplevels(as.data.frame(cc))
  if(nrow(cc)<100L || sum(cc$event_lm)<30L || !all(c("Q1","Q2","Q3","Q4") %in% cc$map_q)) stop("Sparse subgroup")
  coxph(get_formula(covs,exposure),data=cc,ties="efron",x=TRUE)
}
models_all <- list(); warnings_all <- list()
run_model <- function(label,covs,data_list=ddlist,subset_fn=NULL) {
  cat("Model:",label,"\n")
  warnings <- character()
  fits <- tryCatch(withCallingHandlers(lapply(data_list,function(dd){
    if(!is.null(subset_fn)) dd <- dd[subset_fn(dd)]
    fit_one(dd,covs)
  }),warning=function(w){warnings <<- unique(c(warnings,conditionMessage(w)));invokeRestart("muffleWarning")}),error=function(e)e)
  if(inherits(fits,"error")) return(data.table(model=label,term="map_qQ4",status=conditionMessage(fits)))
  if(any(vapply(fits,function(f) any(!is.finite(coef(f))) || any(!is.finite(vcov(f))),logical(1))) || any(grepl("infinite|converge",warnings))) {
    return(data.table(model=label,term="map_qQ4",status=paste("unstable",paste(warnings,collapse=";"))))
  }
  pooled <- if(length(fits)>1) as.data.table(summary(pool(fits),conf.int=TRUE)) else {
    b <- coef(fits[[1]]);s <- sqrt(diag(vcov(fits[[1]])))
    data.table(term=names(b),estimate=b,std.error=s,p.value=2*pnorm(-abs(b/s)),conf.low=b-1.96*s,conf.high=b+1.96*s)
  }
  if("2.5 %" %in% names(pooled) && !"conf.low" %in% names(pooled)) setnames(pooled,c("2.5 %","97.5 %"),c("conf.low","conf.high"))
  pooled[,`:=`(model=label,HR=exp(estimate),lower=exp(conf.low),upper=exp(conf.high),n=fits[[1]]$n,events=fits[[1]]$nevent,m=length(fits),status="estimated")]
  models_all[[label]] <<- fits
  if(length(warnings)) warnings_all[[label]] <<- warnings
  pooled
}
vascular <- c("icd_hypertension","icd_stroke")
context_vars <- c("careunit_group","admission_group","cardio_proc_device")
phenotype <- c(primary,vascular,context_vars)
extended <- c(phenotype,"advanced_monitor","ne_log1p","service_group")
model_specs <- list(
 "Primary F52 reproduction"=primary,
 "F52 phenotype reproduction"=phenotype,
 "Primary plus vascular history"=c(primary,vascular),
 "Primary plus ICU type"=c(primary,"careunit_group"),
 "Primary plus admission type"=c(primary,"admission_group"),
 "Primary plus procedures and devices"=c(primary,"cardio_proc_device","advanced_monitor"),
 "Primary plus NE dose"=c(primary,"ne_log1p"),
 "Primary plus hospital service"=c(primary,"service_group"),
 "Extended clinical context"=extended,
 "SOFA domain replacement"=c(setdiff(extended,"sofa"),domains),
 "Noncardiovascular SOFA"=c(setdiff(extended,"sofa"),"sofa_noncv"),
 "Non-CNS SOFA"=c(setdiff(extended,"sofa"),"sofa_noncns"),
 "Extended without MV"=setdiff(extended,"itvtn_24h_vent_tag"),
 "Extended without RRT"=setdiff(extended,"itvtn_24h_rrt_tag"),
 "Extended without sedation"=setdiff(extended,"drug_24h_sedative_tag")
)
needed <- unique(unlist(model_specs))
fwrite(data.table(variable=needed,missing=vapply(needed,function(v)sum(is.na(ddlist[[1]][[v]])),integer(1))),"04_QC/covariate_missingness.csv")
common <- ddlist[[1]][complete.cases(ddlist[[1]][,needed,with=FALSE]),stay_id]
common_list <- lapply(ddlist,function(d)d[stay_id %in% common])
res <- rbindlist(lapply(names(model_specs),function(nm)run_model(nm,model_specs[[nm]],common_list)),fill=TRUE)
res[,sample_scope:="Observed six-domain SOFA complete cases"]
fwrite(res,"02_Results/P01_block_models_all_coefficients.csv")
q4 <- res[term=="map_qQ4"]
reference <- q4[model=="Primary F52 reproduction",estimate]
q4[,descriptive_logHR_change_pct:=100*(reference-estimate)/reference]
fwrite(q4,"02_Results/P02_Q4_coefficient_changes.csv")
reproduce <- rbindlist(list(run_model("Primary all eligible",primary),run_model("Phenotype all eligible",phenotype)),fill=TRUE)
fwrite(reproduce,"02_Results/P00_F52_reproduction.csv")

strata <- list("All"=function(d)rep(TRUE,nrow(d)),"Medical ICU"=function(d)d$careunit_group=="Medical ICU",
 "Cardiac ICU"=function(d)d$careunit_group=="Cardiac ICU","Surgical/Trauma ICU"=function(d)d$careunit_group=="Surgical/Trauma ICU",
 "Neuro ICU/Stepdown"=function(d)d$careunit_group=="Neuro ICU/Stepdown",
 "Exclude neuro ICU"=function(d)d$careunit_group!="Neuro ICU/Stepdown",
 "Exclude surgical/trauma ICU"=function(d)d$careunit_group!="Surgical/Trauma ICU",
 "Exclude cardio procedures/devices"=function(d)d$cardio_proc_device==0,
 "Sepsis established by 24 h"=function(d)d$sepsis_by24,
 "No vasopressor"=function(d)d$ne_equiv_mean_0_24h==0,
 "Vasopressor"=function(d)d$ne_equiv_mean_0_24h>0,
 "No MV"=function(d)d$itvtn_24h_vent_tag==0,"MV"=function(d)d$itvtn_24h_vent_tag==1,
 "No RRT"=function(d)d$itvtn_24h_rrt_tag==0,"RRT"=function(d)d$itvtn_24h_rrt_tag==1)
strat <- rbindlist(lapply(names(strata),function(nm)run_model(nm,primary,subset_fn=strata[[nm]])),fill=TRUE)
fwrite(strat,"02_Results/V02_MI_subgroup_models.csv")
counts <- rbindlist(lapply(names(strata),function(nm){d<-ddlist[[1]][strata[[nm]](ddlist[[1]])];d[,.(n=.N,events=sum(event_lm)),by=.(map_q)][,stratum:=nm]}))
fwrite(counts,"02_Results/V03_subgroup_event_counts.csv")

# Temporal stability is based on observed hours, with no cross-period imputation.
obs <- copy(dat)
obs[,map_q:=relevel(cut(mean_observed,c(-Inf,71.40,75.72,80.99,Inf),labels=c("Q1","Q2","Q3","Q4"),right=FALSE),"Q2")]
obs[,gender:=factor(gender)]
temporal <- list()
for(per in c("Earlier <=2016","Later >=2017","Boundary interval")) {
 for(covname in c("Primary","Phenotype")) {
  temp <- obs[period==per]
  temporal[[paste(per,covname)]] <- run_model(paste(per,covname),if(covname=="Primary")primary else phenotype,list(temp))
 }
 temporal[[paste(per,"complete24")]] <- run_model(paste(per,"complete24"),primary,list(obs[period==per & n_hour_obs==24]))
}
fwrite(rbindlist(temporal,fill=TRUE),"02_Results/V04_temporal_observed_models.csv")
service <- rbindlist(lapply(levels(dat$service_group),function(s)run_model(s,primary,subset_fn=function(d)d$service_group==s)),fill=TRUE)
fwrite(service,"02_Results/V05_service_models.csv")

# Full-cohort interactions, descriptions and PH diagnostics are in stage 07.
if(length(warnings_all)) capture.output(warnings_all,file="04_QC/model_warnings.txt")
capture.output(sessionInfo(),file="04_QC/R_analysis_sessionInfo.txt")
cat("Finished",format(Sys.time()),"\n")
sink();close(logfile)
