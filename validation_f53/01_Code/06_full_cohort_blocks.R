invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(data.table);library(mice);library(survival)})
spec<-readRDS("restricted_cache/model_spec.rds")
primary<-spec$covars
domains<-paste0("sofa_",c("respiration","coagulation","liver","cardiovascular","cns","renal"))
pheno<-c(primary,"icd_hypertension","icd_stroke","careunit_group","admission_group","cardio_proc_device")
extended<-c(pheno,"advanced_monitor","ne_log1p","service_group")
keep<-unique(c("stay_id","time_lm_days","event_lm","map_q",extended,domains,"sepsis_by24","n_hour_obs","mean_map",spec$hour_cols))
ddlist<-lapply(readRDS("restricted_cache/mimic_m40_enriched.rds"),function(d)d[,..keep])
missing<-data.table(domain=domains,n_missing=vapply(domains,function(v)sum(is.na(ddlist[[1]][[v]])),integer(1)),n=nrow(ddlist[[1]]))
fwrite(missing,"04_QC/SOFA_raw_component_missingness.csv")
for(i in seq_along(ddlist)) {
 d<-ddlist[[i]]
 d[,domains_complete:=complete.cases(.SD),.SDcols=domains]
 # Match the published first_day_sofa total-score convention. Preserve missingness above.
 for(v in domains) set(d,which(is.na(d[[v]])),v,0)
 stopifnot(all(d$sofa==rowSums(d[,..domains])))
 d[,sofa_noncv:=sofa-sofa_cardiovascular]
 d[,sofa_noncns:=sofa-sofa_cns]
 ddlist[[i]]<-d
}
models_all<-list();warnings_all<-list()
expr<-as.list(parse("01_Code/03_phenotype_validation.R"))
for(nm in c("get_formula","fit_one","run_model")) {
 x<-Filter(function(x)is.call(x)&&identical(x[[1]],as.name("<-"))&&identical(x[[2]],as.name(nm)),expr)
 stopifnot(length(x)==1);eval(x[[1]])
}
specs<-list("Primary model"=primary,"Plus vascular history"=c(primary,"icd_hypertension","icd_stroke"),
 "Plus ICU type"=c(primary,"careunit_group"),"Plus admission type"=c(primary,"admission_group"),
 "Plus procedures and devices"=c(primary,"cardio_proc_device","advanced_monitor"),"Plus NE dose"=c(primary,"ne_log1p"),
 "Plus hospital service"=c(primary,"service_group"),"F52 phenotype model"=pheno,
 "Extended clinical context"=extended,"SOFA domain replacement"=c(setdiff(extended,"sofa"),domains),
 "Noncardiovascular SOFA"=c(setdiff(extended,"sofa"),"sofa_noncv"),"Non-CNS SOFA"=c(setdiff(extended,"sofa"),"sofa_noncns"),
 "Extended without MV"=setdiff(extended,"itvtn_24h_vent_tag"),"Extended without RRT"=setdiff(extended,"itvtn_24h_rrt_tag"),
 "Extended without sedation"=setdiff(extended,"drug_24h_sedative_tag"))
vars<-unique(unlist(specs))
ids<-ddlist[[1]][complete.cases(ddlist[[1]][,vars,with=FALSE]),stay_id]
ddlist<-lapply(ddlist,function(d)d[stay_id %in% ids])
out<-rbindlist(lapply(names(specs),function(nm)run_model(nm,specs[[nm]])),fill=TRUE)
out[,sample_scope:="Common full cohort; SOFA missing components scored zero per source total"]
fwrite(out,"02_Results/P07_full_cohort_models.csv")
q<-out[term=="map_qQ4"]
ref<-q[model=="Primary model",estimate]
q[,descriptive_logHR_change_pct:=100*(ref-estimate)/ref]
fwrite(q,"02_Results/P08_full_cohort_Q4_changes.csv")
early<-rbindlist(list(run_model("Early sepsis primary",primary,subset_fn=function(d)d$sepsis_by24),
 run_model("Early sepsis extended context",extended,subset_fn=function(d)d$sepsis_by24),
 run_model("Complete24 primary",primary,subset_fn=function(d)d$n_hour_obs==24),
 run_model("Complete24 extended context",extended,subset_fn=function(d)d$n_hour_obs==24)),fill=TRUE)
fwrite(early,"02_Results/P09_early_sepsis_complete24_context.csv")
saveRDS(ddlist,"restricted_cache/mimic_m40_compact_scored_domains.rds")
saveRDS(models_all,"restricted_cache/full_cohort_fits.rds")
cat("Full cohort explanation complete\n")
