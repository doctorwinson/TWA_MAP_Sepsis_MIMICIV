invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(data.table);library(mice);library(survival)})
spec<-readRDS("restricted_cache/model_spec.rds");primary<-spec$covars
dl<-readRDS("restricted_cache/mimic_m40_compact_scored_domains.rds")
for(i in seq_along(dl))dl[[i]][,vaso_any:=as.integer(ne_log1p>0)]
ints<-list()
for(v in c("careunit_group","service_group","itvtn_24h_vent_tag","itvtn_24h_rrt_tag","vaso_any")) {
  groups<-dl[[1]][,.(n=.N,events=sum(event_lm)),by=v]
  valid<-groups[n>=50 & events>=10,get(v)]
  subset_dl<-lapply(dl,function(d)droplevels(d[get(v) %in% valid]))
  ws<-character()
  row<-tryCatch(withCallingHandlers({
    rhs<-unique(c("map_q",primary,v));f0<-as.formula(paste("Surv(time_lm_days,event_lm) ~",paste(rhs,collapse=" + ")))
    f1<-update(f0,paste(". ~ . + map_q:",v))
    fits0<-lapply(subset_dl,function(d)coxph(f0,data=d))
    fits1<-lapply(subset_dl,function(d)coxph(f1,data=d))
    if(any(vapply(fits1,function(f)anyNA(coef(f)),logical(1))))stop("Aliased coefficient")
    if(any(grepl("infinite|converge",ws)))stop("Unstable interaction fit")
    z<-as.matrix(mice::D1(as.mira(fits1),as.mira(fits0))$result)
    pc<-grep("^P\\(>|p.value",colnames(z),ignore.case=TRUE)
    stopifnot(length(pc)==1L)
    data.table(variable=v,p.value=as.numeric(z[1,pc]),n=nrow(subset_dl[[1]]),status="estimated; levels with <50 patients or <10 events excluded")
  },warning=function(w){ws<<-c(ws,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)data.table(variable=v,p.value=NA_real_,n=nrow(subset_dl[[1]]),status=conditionMessage(e)))
  ints[[v]]<-row
}
d<-as.data.table(readRDS("restricted_cache/mimic_enriched_preimputation.rds"))
d[,map_q:=relevel(cut(mean_observed,c(-Inf,71.4,75.72,80.99,Inf),labels=c("Q1","Q2","Q3","Q4"),right=FALSE),"Q2")]
d<-droplevels(d[period!="Boundary interval"])
d<-d[complete.cases(d[,c(primary,"event_lm","time_lm_days","map_q","period"),with=FALSE])]
f0<-as.formula(paste("Surv(time_lm_days,event_lm) ~ map_q + period +",paste(primary,collapse=" + ")))
f1<-update(f0,.~.+map_q:period)
a<-anova(coxph(f0,data=d),coxph(f1,data=d),test="Chisq")
ints[["period"]]<-data.table(variable="period_observed_hours",p.value=a[2,"Pr(>|Chi|)"],n=nrow(d),status="single-data likelihood ratio interaction; boundary intervals excluded")
it<-rbindlist(ints,fill=TRUE);it[,p_BH:=p.adjust(p.value,"BH")]
fwrite(it,"02_Results/V06_pooled_interactions.csv")

d1<-dl[[1]]
catvars<-c("careunit_group","service_group","admission_group","cardio_proc_device","advanced_monitor","itvtn_24h_vent_tag","itvtn_24h_rrt_tag","drug_24h_sedative_tag","vaso_any")
desc<-rbindlist(lapply(catvars,function(v){x<-d1[,.N,by=c("map_q",v)];setnames(x,c(v,"N"),c("level","n"));x[,level:=as.character(level)];x[,denominator:=sum(n),by=map_q];x[,`:=`(variable=v,percent=100*n/denominator)];x}),fill=TRUE)
fwrite(desc,"02_Results/P03_phenotype_distribution.csv")
dom<-paste0("sofa_",c("respiration","coagulation","liver","cardiovascular","cns","renal"))
descnum<-rbindlist(lapply(c(dom,"ne_log1p","sofa"),function(v)d1[,.(n=sum(!is.na(get(v))),median=as.numeric(median(get(v),na.rm=TRUE)),q25=as.numeric(quantile(get(v),.25,na.rm=TRUE)),q75=as.numeric(quantile(get(v),.75,na.rm=TRUE))),by=map_q][,variable:=v]))
fwrite(descnum,"02_Results/P04_SOFA_NE_distribution.csv")
fits<-readRDS("restricted_cache/full_cohort_fits.rds")
ph<-rbindlist(lapply(c("Primary model","F52 phenotype model","Extended clinical context"),function(nm){z<-cox.zph(fits[[nm]][[1]],terms=TRUE)$table;data.table(model=nm,term=rownames(z),chisq=z[,1],df=z[,2],p.value=z[,3],scope="Imputation 1 diagnostic; not a pooled test")}))
fwrite(ph,"04_QC/proportional_hazards_imputation1.csv")

# Separate early/late follow-up associations address a possible PH violation.
tfits<-lapply(dl,function(d){
 z<-survSplit(Surv(time_lm_days,event_lm)~.,data=as.data.frame(d[,c("time_lm_days","event_lm","map_q",primary),with=FALSE]),cut=6,episode="window")
 z$late<-as.integer(z$window==2)
 for(q in c("Q1","Q3","Q4"))z[[paste0(q,"_late")]]<-as.integer(z$map_q==q)*z$late
 f<-as.formula(paste("Surv(tstart,time_lm_days,event_lm) ~ map_q + Q1_late + Q3_late + Q4_late + strata(window) +",paste(primary,collapse=" + ")))
 coxph(f,data=z)
})
tr<-as.data.table(summary(pool(tfits),conf.int=TRUE));fwrite(tr,"02_Results/V07_time_interaction_coefficients.csv")
contrasts<-rbindlist(lapply(c("Days 1-7 after ICU admission","Days 7-30 after ICU admission"),function(label){
 b<-vapply(tfits,function(f){x<-coef(f);x["map_qQ4"]+if(grepl("7-30",label))x["Q4_late"]else 0},numeric(1))
 u<-vapply(tfits,function(f){v<-vcov(f);if(grepl("7-30",label))v["map_qQ4","map_qQ4"]+v["Q4_late","Q4_late"]+2*v["map_qQ4","Q4_late"]else v["map_qQ4","map_qQ4"]},numeric(1))
 z<-pool.scalar(b,u,n=nrow(dl[[1]]),k=length(coef(tfits[[1]])))
 data.table(window=label,HR=exp(z$qbar),lower=exp(z$qbar-qt(.975,z$df)*sqrt(z$t)),upper=exp(z$qbar+qt(.975,z$df)*sqrt(z$t)),m=40L)
}))
fwrite(contrasts,"02_Results/V08_time_specific_Q4.csv")

mi<-readRDS("restricted_cache/mice40_full.rds")$imp
cm<-as.data.table(as.table(mi$chainMean));setnames(cm,c("hour","iteration","chain","mean"))
cm[,`:=`(iteration=as.integer(as.character(iteration)),chain=as.integer(gsub("[^0-9]","",as.character(chain))))]
fwrite(cm,"04_QC/mice_chain_means.csv")
diags<-rbindlist(lapply(spec$hour_cols,function(h){
 observed<-mi$data[[h]][!is.na(mi$data[[h]])];imputed<-unlist(mi$imp[[h]])
 data.table(hour=h,n_observed=length(observed),n_missing=sum(is.na(mi$data[[h]])),obs_q01=quantile(observed,.01),obs_median=median(observed),obs_q99=quantile(observed,.99),imp_q01=quantile(imputed,.01),imp_median=median(imputed),imp_q99=quantile(imputed,.99),imp_min=min(imputed),imp_max=max(imputed))
}))
fwrite(diags,"04_QC/imputed_value_distribution.csv")
capture.output(sessionInfo(),file="04_QC/R_diagnostics_sessionInfo.txt")
cat("Diagnostics and interactions complete\n")
