invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(data.table);library(survival);library(mice)})
sp<-readRDS("restricted_cache/model_spec.rds");a<-readRDS("restricted_cache/analysis_specs.rds")
w<-as.data.table(readRDS("restricted_cache/weight_provenance.rds"))
d<-merge(copy(a$raw),w[,.(stay_id,official_day_weight,official_age,age_score)],by="stay_id",all.x=TRUE)
d[,map_q:=relevel(cut(mean_observed,c(-Inf,sp$cutpoints[2:4],Inf),labels=c("Q1","Q2","Q3","Q4"),right=TRUE),"Q2")]
d<-droplevels(d[!is.na(official_day_weight)])
out<-rbindlist(lapply(c("Source weight, paired sample","Official first-day weight, paired sample"),function(nm){
 dd<-copy(d);if(grepl("Official",nm))dd[,weight:=official_day_weight]
 f<-coxph(as.formula(paste("Surv(time_lm_days,event_lm)~map_q+",paste(a$primary,collapse="+"))),data=dd)
 b<-coef(f);se<-sqrt(diag(vcov(f)))
 data.table(model=nm,term=names(b),HR=exp(b),lower=exp(b-1.96*se),upper=exp(b+1.96*se),p.value=2*pnorm(-abs(b/se)),n=f$n,events=f$nevent,exposure="Observed-hour MAP; no imputation")
}))
fwrite(out,"02_Results/weight_source_sensitivity.csv")
x<-as.data.table(readRDS("restricted_cache/sepsis_source_verification.rds"))
x<-x[stay_id %in% a$raw$stay_id]
x[,death_exact:=deathtime]
x[is.na(death_exact),death_exact:=as.POSIXct(dod,tz="UTC")]
x[,death_hours:=as.numeric(difftime(death_exact,intime,units="hours"))]
x[,event_direct:=as.integer(!is.na(death_hours)&death_hours>24&death_hours<=720)]
fwrite(data.table(metric=c("rows","event_disagreement_vs_timestamp_or_date","event_disagreement_among_exact_hospital_deaths","date_only_deaths"),value=c(nrow(x),sum(x$event_direct!=x$day30_outcome,na.rm=TRUE),sum(x$event_direct!=x$day30_outcome & !is.na(x$deathtime),na.rm=TRUE),sum(is.na(x$deathtime)&!is.na(x$dod)))),"04_QC/direct_outcome_audit.csv")
print(fread("04_QC/direct_outcome_audit.csv"))
