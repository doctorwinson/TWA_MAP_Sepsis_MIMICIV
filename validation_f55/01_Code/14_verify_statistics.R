invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(data.table);library(mice);library(survival)})
checks<-list()
check<-function(label,ok){checks[[length(checks)+1L]]<<-data.table(check=label,pass=isTRUE(ok));if(!isTRUE(ok))cat("FAIL",label,"\n")}
src<-as.data.table(readRDS("restricted_cache/source_enriched.rds"))
eligible<-as.data.table(readRDS("restricted_cache/eligible_primary.rds"))
sp<-readRDS("restricted_cache/model_spec.rds")
dl<-readRDS("restricted_cache/analysis_m40.rds")
N<-nrow(eligible);E<-sum(eligible$event_lm)
check("source contains analysis cohort",all(eligible$stay_id %in% src$stay_id))
check("all analysis stays have onset by24",all(eligible$onset_hours<=24))
check("unique patients",!anyDuplicated(eligible$subject_id))
check("all observed hours >=18",all(eligible$n_hour_obs>=18))
check("all at-risk follow-up positive and <=29 days",all(eligible$time_lm_days>0 & eligible$time_lm_days<=29))
check("all ICU occupancy >=24h",all(eligible$icu_hours>=24))
for(scope in c("primary","enriched")){
 imp<-readRDS(paste0("restricted_cache/mice40_",scope,".rds"))
 check(paste(scope,"m40 iteration30"),imp$m==40 && imp$iteration==30)
 check(paste(scope,"sex factor"),is.factor(imp$data$gender))
 check(paste(scope,"sex used as MAP predictor"),all(imp$predictorMatrix[sp$hours,"gender"]==1))
 for(i in 1:40){
  d<-complete(imp,i)
  check(paste(scope,i,"complete and physiologic MAP"),!anyNA(d[,sp$hours]) && min(as.matrix(d[,sp$hours]))>=20 && max(as.matrix(d[,sp$hours]))<=200)
  same<-vapply(sp$hours,function(h){ix<-!is.na(imp$data[[h]]);identical(d[[h]][ix],imp$data[[h]][ix])},logical(1))
  check(paste(scope,i,"observed unchanged"),all(same))
 }
 conv<-as.data.table(mice::convergence(imp))
 fwrite(conv,paste0("04_QC/MI_",scope,"_convergence_final.csv"))
 h<-conv[.it==30 & vrb %in% sp$hours]
 fwrite(h,paste0("04_QC/MI_",scope,"_convergence_last.csv"))
 check(paste(scope,"finite hour chain diagnostics"),all(is.finite(h$psrf)))
}
check("40 completed analysis datasets",length(dl)==40)
check("events unchanged across imputations",all(vapply(dl,function(d)sum(d$event_lm)==E,logical(1))))
for(d in dl)check("fixed quartile categorization",identical(as.character(d$map_q),as.character(cut(d$mean_map,c(-Inf,sp$cutpoints[2:4],Inf),labels=c("Q1","Q2","Q3","Q4"),right=TRUE))))
counts<-fread("02_Results/mortality_by_quartile.csv")
check("descriptive count total",sum(counts$n)==N && sum(counts$events)==E)
truth<-dl[[1]][,.(n=.N,events=sum(event_lm)),by=map_q]
c<-merge(counts,truth,by="map_q")
check("quartile counts match first completed dataset",all(c$n.x==c$n.y & c$events.x==c$events.y))
models<-fread("02_Results/main_models.csv")
fits<-readRDS("restricted_cache/core_fits.rds")
for(nm in names(fits)){
 p<-as.data.table(summary(pool(fits[[nm]]),conf.int=TRUE,exponentiate=TRUE))
 x<-merge(models[model==nm],p,by="term")
 check(paste(nm,"independent pool point estimates"),max(abs(x$HR-x$estimate.y))<1e-10)
 check(paste(nm,"common sample size"),all(x$n==N & x$events==E))
}
for(file in c("main_models.csv","sensitivity_models.csv","subgroup_models.csv","temporal_models.csv","enriched_MI_models.csv")){
 d<-fread(file.path("02_Results",file))
 if("status" %in% names(d)){
  failed<-d[status!="estimated"]
  if(nrow(failed))fwrite(failed,paste0("04_QC/nonestimable_",file))
  check(paste(file,"no unreported model failure"),nrow(failed)==0 || (file=="temporal_models.csv" && nrow(failed)==1 && grepl("Boundary",failed$model,ignore.case=TRUE) && grepl("Extended",failed$model) && grepl("infinite",failed$status) && failed$events==27))
  d<-d[status=="estimated"]
 }
 check(paste(file,"ordered positive finite confidence limits"),all(is.finite(d$HR)&is.finite(d$lower)&is.finite(d$upper)&d$lower>0&d$lower<=d$HR&d$upper>=d$HR))
}
curves<-fread("02_Results/spline_curves.csv")
check("spline reference HR=1",all(curves[abs(x-reference)<1e-12,HR]==1))
check("spline limits ordered",all(curves$lower<=curves$HR & curves$upper>=curves$HR))
check("interaction BH not below raw P",all(fread("02_Results/interactions.csv")[,p_BH>=p.value]))
flow<-fread("02_Results/F55_flow.csv")
check("flow end before covariate exclusion",flow[stage=="At least18 observed hours",n]==nrow(src))
check("flow final denominator",tail(flow$n,1)==N)
out<-rbindlist(checks);fwrite(out,"04_QC/numerical_verification.csv")
cat(sum(out$pass),"/",nrow(out),"checks passed\n")
if(!all(out$pass))stop("Numerical checks failed")
