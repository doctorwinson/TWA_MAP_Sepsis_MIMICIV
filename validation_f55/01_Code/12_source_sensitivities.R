invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(data.table);library(survival);library(mice)})
sp<-readRDS("restricted_cache/model_spec.rds");a<-readRDS("restricted_cache/analysis_specs.rds")
src<-as.data.table(readRDS("restricted_cache/source_enriched.rds"))
d<-merge(copy(a$raw),src[,.(stay_id,weight_admit,invasive_or_trach)],by="stay_id",all.x=TRUE)
d[,map_q:=relevel(cut(mean_observed,c(-Inf,sp$cutpoints[2:4],Inf),labels=c("Q1","Q2","Q3","Q4"),right=TRUE),"Q2")]
d<-droplevels(d[!is.na(weight_admit)])
out<-rbindlist(lapply(c("First-day mean weight, paired sample","Admission weight, paired sample","Invasive ventilation or tracheostomy, paired sample"),function(nm){
 dd<-copy(d)
 if(nm=="Admission weight, paired sample")dd[,weight:=weight_admit]
 if(grepl("tracheostomy",nm))dd[,itvtn_24h_vent_tag:=invasive_or_trach]
 f<-coxph(as.formula(paste("Surv(time_lm_days,event_lm)~map_q+",paste(a$primary,collapse="+"))),data=dd)
 b<-coef(f);se<-sqrt(diag(vcov(f)))
 data.table(model=nm,term=names(b),HR=exp(b),lower=exp(b-1.96*se),upper=exp(b+1.96*se),p.value=2*pnorm(-abs(b/se)),n=f$n,events=f$nevent,exposure="Observed-hour MAP; no imputation")
}))
fwrite(out,"02_Results/weight_source_sensitivity.csv")
