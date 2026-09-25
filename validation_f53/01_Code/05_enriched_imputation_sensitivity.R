invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(data.table);library(mice);library(survival)})
source("01_Code/legacy_mi_function.R")
full_scope<-identical(Sys.getenv("F53_IMPUTATION_SCOPE","complete_domains"),"full")
tag<-if(full_scope)"enriched_full"else "enriched_complete_domains"
logfile<-file(paste0("04_QC/",tag,".log"),"wt");sink(logfile)
dat<-as.data.table(readRDS("restricted_cache/mimic_enriched_preimputation.rds"))
spec<-readRDS("restricted_cache/model_spec.rds")
primary<-spec$covars
pheno<-c(primary,"icd_hypertension","icd_stroke","careunit_group","admission_group","cardio_proc_device")
extended<-c(pheno,"advanced_monitor","ne_log1p","service_group")
domains<-paste0("sofa_",c("respiration","coagulation","liver","cardiovascular","cns","renal"))
predictors<-setdiff(unique(c(extended,domains,"sepsis_by24")),"sofa")
if(full_scope) {
 for(v in domains) {
  mv<-paste0(v,"_missing");dat[,(mv):=as.integer(is.na(get(v)))];predictors<-c(predictors,mv)
  set(dat,which(is.na(dat[[v]])),v,0)
 }
 stopifnot(all(dat$sofa==rowSums(dat[,..domains])))
}
keep<-complete.cases(dat[,predictors,with=FALSE])
dat<-dat[keep]
dat[,gender:=factor(gender)]
cat("Fixed complete-predictor cohort",nrow(dat),"\n")
res<-run_mi_hours_only(dat,predictors,spec$hour_cols,m=40,maxit=10,seed=42)
saveRDS(res,paste0("restricted_cache/mice40_",tag,".rds"))
data_list<-lapply(1:40,function(i){
 d<-as.data.table(complete(res$imp,i));d<-merge(d,dat[,.(stay_id,sofa)],by="stay_id",all.x=TRUE)
 d[,mean_map:=rowMeans(.SD),.SDcols=spec$hour_cols]
 d[,map_q:=relevel(cut(mean_map,c(-Inf,spec$cutpoints[2:4],Inf),labels=c("Q1","Q2","Q3","Q4")),"Q2")];d
})
specs<-list("Enriched MI primary"=primary,"Enriched MI phenotype"=pheno,"Enriched MI extended"=extended,
 "Enriched MI SOFA components"=c(setdiff(extended,"sofa"),domains),"Enriched MI primary early sepsis"=primary)
out<-rbindlist(lapply(names(specs),function(nm){
 fits<-lapply(data_list,function(d){
   if(nm=="Enriched MI primary early sepsis")d<-d[sepsis_by24==TRUE]
   coxph(as.formula(paste("Surv(time_lm_days,event_lm) ~ map_q +",paste(specs[[nm]],collapse=" + "))),data=d)
 })
 z<-as.data.table(summary(pool(fits),conf.int=TRUE,exponentiate=TRUE))
 z[,`:=`(model=nm,n=fits[[1]]$n,events=fits[[1]]$nevent,m=40L)];z
}),fill=TRUE)
out[,scope:=if(full_scope)"Full cohort; scored-zero SOFA domains with missingness indicators"else "Observed six-domain SOFA complete cases"]
fwrite(out,paste0("02_Results/",if(full_scope)"P10"else "P05","_",tag,"_MI_sensitivity.csv"))
if(!is.null(res$imp$loggedEvents))fwrite(as.data.table(res$imp$loggedEvents),paste0("04_QC/",tag,"_logged_events.csv"))
cat("Finished",format(Sys.time()),"\n");sink();close(logfile)
cat("Phenotype-enriched m=40 sensitivity completed\n")
