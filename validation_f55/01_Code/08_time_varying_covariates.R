invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(data.table);library(mice);library(survival)})
dl<-readRDS("restricted_cache/analysis_m40.rds");sp<-readRDS("restricted_cache/analysis_specs.rds")
fs<-lapply(dl,function(d){
 mm<-model.matrix(as.formula(paste("~map_q +",paste(sp$primary,collapse=" + "))),data=d)[,-1,drop=FALSE]
 z<-data.frame(time=d$time_lm_days,event=d$event_lm,mm,check.names=FALSE)
 z<-survSplit(Surv(time,event)~.,data=z,cut=6,episode="window")
 cn<-colnames(mm)
 for(v in cn)z[[paste0(v,"_late")]]<-z[[v]]*as.integer(z$window==2)
 f<-as.formula(paste("Surv(tstart,time,event) ~",paste(c(cn,paste0(cn,"_late"),"strata(window)"),collapse=" + ")))
 fit<-coxph(f,data=z);stopifnot(all(is.finite(coef(fit))));fit
})
result<-rbindlist(lapply(c(FALSE,TRUE),function(late){
 b<-vapply(fs,function(f)coef(f)["map_qQ4"]+if(late)coef(f)["map_qQ4_late"]else 0,numeric(1))
 u<-vapply(fs,function(f){v<-vcov(f);v["map_qQ4","map_qQ4"]+if(late)v["map_qQ4_late","map_qQ4_late"]+2*v["map_qQ4","map_qQ4_late"]else 0},numeric(1))
 z<-pool.scalar(b,u,n=fs[[1]]$nevent,k=length(coef(fs[[1]])))
 data.table(window=if(late)"ICU days 7-30"else"ICU days 1-7",HR=exp(z$qbar),lower=exp(z$qbar-qt(.975,z$df)*sqrt(z$t)),upper=exp(z$qbar+qt(.975,z$df)*sqrt(z$t)),adjustment="All primary covariate and exposure coefficients allowed to differ by window")
}))
fwrite(result,"02_Results/time_varying_covariates_Q4.csv")
print(result)
