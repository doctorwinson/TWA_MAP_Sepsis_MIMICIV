invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(data.table);library(mice);library(survival);library(splines)})
sp<-readRDS("restricted_cache/model_spec.rds");asp<-readRDS("restricted_cache/analysis_specs.rds")
dl<-readRDS("restricted_cache/analysis_m40.rds");primary<-asp$primary
fits<-readRDS("restricted_cache/core_fits.rds")
ph<-rbindlist(lapply(names(fits),function(nm)rbindlist(lapply(seq_along(fits[[nm]]),function(i){
 z<-cox.zph(fits[[nm]][[i]],terms=TRUE)$table
 data.table(model=nm,imputation=i,term=rownames(z),chisq=z[,1],df=z[,2],p.value=z[,3])
}))))
fwrite(ph,"02_Results/PH_all_imputations.csv")
fwrite(ph[,.(min_p=min(p.value),median_p=median(p.value),max_p=max(p.value),n_p_below05=sum(p.value<.05)),by=.(model,term)],"02_Results/PH_summary.csv")

test_p<-function(x){z<-as.matrix(x$result);j<-grep("^P\\(>|p.value",colnames(z),ignore.case=TRUE);stopifnot(length(j)==1);as.numeric(z[1,j])}
splineres<-list();splinetests<-list()
for(exposure in c("mean_map","ttr65")) {
 ref<-if(exposure=="mean_map")sp$reference else median(dl[[1]]$ttr65)
 knots<-if(exposure=="mean_map")sp$knots else c(1/24,4/24,12/24)
 boundary<-if(exposure=="mean_map")sp$boundary else c(0,1)
 basis<-function(x)ns(x,knots=knots,Boundary.knots=boundary)
 refb<-as.numeric(basis(ref))
 g<-seq(boundary[1],boundary[2],length.out=101)
 centered<-sweep(basis(g),2,refb,"-")
 qr0<-qr(cbind(g-ref,centered),tol=1e-8)
 stopifnot(qr0$rank==4,qr0$pivot[1]==1)
 take<-qr0$pivot[2:4]-1
 addbasis<-function(d){d<-copy(d);bb<-basis(d[[exposure]]);for(j in 1:3)d[,(paste0("nl",j)):=bb[,take[j]]];d}
 dd<-lapply(dl,addbasis)
 for(adj in c("Primary","Extended")){
  covs<-if(adj=="Primary")primary else asp$extended
  f0<-as.formula(paste("Surv(time_lm_days,event_lm) ~",paste(c(covs,exposure),collapse=" + ")))
  f1<-update(f0,.~.+nl1+nl2+nl3)
  f2<-as.formula(paste("Surv(time_lm_days,event_lm) ~",paste(covs,collapse=" + ")))
  fs<-lapply(dd,function(d)coxph(f1,data=d))
  fl<-lapply(dd,function(d)coxph(f0,data=d))
  fn<-lapply(dd,function(d)coxph(f2,data=d))
  splinetests[[paste(exposure,adj)]]<-data.table(exposure=exposure,adjustment=adj,nonlinear_p=test_p(D1(as.mira(fs),as.mira(fl))),overall_p=test_p(D1(as.mira(fs),as.mira(fn))),reference=ref,knots=paste(knots,collapse=";"),n=nrow(dd[[1]]))
  xr<-as.numeric(quantile(dl[[1]][[exposure]],c(.01,.99)))
  xx<-sort(unique(c(seq(xr[1],xr[2],length.out=201),ref)))
  bb<-sweep(basis(xx),2,refb,"-")[,take,drop=FALSE]
  contrasts<-cbind(xx-ref,bb);colnames(contrasts)<-c(exposure,"nl1","nl2","nl3")
  rr<-rbindlist(lapply(seq_along(xx),function(i){
   a<-contrasts[i,];nm<-names(a)
   if(all(abs(a)<1e-12))return(data.table(x=xx[i],HR=1,lower=1,upper=1))
   b<-vapply(fs,function(f)sum(a*coef(f)[nm]),numeric(1))
   v<-vapply(fs,function(f)as.numeric(t(a)%*%vcov(f)[nm,nm]%*%a),numeric(1))
   z<-pool.scalar(b,v,n=fs[[1]]$nevent,k=length(coef(fs[[1]])))
   data.table(x=xx[i],HR=exp(z$qbar),lower=exp(z$qbar-qt(.975,z$df)*sqrt(z$t)),upper=exp(z$qbar+qt(.975,z$df)*sqrt(z$t)))
  }))
  rr[,`:=`(exposure=exposure,adjustment=adj,reference=ref)]
  splineres[[paste(exposure,adj)]]<-rr
 }
}
fwrite(rbindlist(splineres),"02_Results/spline_curves.csv")
fwrite(rbindlist(splinetests),"02_Results/spline_tests.csv")

ints<-list()
for(v in c("careunit_group","service_group","itvtn_24h_vent_tag","itvtn_24h_rrt_tag","vaso_any")) {
 groups<-dl[[1]][,.(n=.N,events=sum(event_lm)),by=v]
 valid<-groups[n>=50 & events>=10,get(v)]
 dd<-lapply(dl,function(d)droplevels(d[get(v) %in% valid]))
 row<-tryCatch({
  rhs<-unique(c("map_q",primary,v));f0<-as.formula(paste("Surv(time_lm_days,event_lm) ~",paste(rhs,collapse=" + ")))
  f1<-update(f0,paste(". ~ . + map_q:",v))
  f0s<-lapply(dd,function(d)coxph(f0,data=d));f1s<-lapply(dd,function(d)coxph(f1,data=d))
  stopifnot(all(vapply(f1s,function(f)all(is.finite(coef(f))),logical(1))))
  data.table(variable=v,p.value=test_p(D1(as.mira(f1s),as.mira(f0s))),n=nrow(dd[[1]]),status="estimated")
 },error=function(e)data.table(variable=v,p.value=NA_real_,n=nrow(dd[[1]]),status=conditionMessage(e)))
 ints[[v]]<-row
}
obs<-copy(asp$raw)[period!="Boundary interval"]
obs[,map_q:=relevel(cut(mean_observed,c(-Inf,sp$cutpoints[2:4],Inf),labels=c("Q1","Q2","Q3","Q4"),right=TRUE),"Q2")]
f0<-as.formula(paste("Surv(time_lm_days,event_lm) ~ map_q + period +",paste(primary,collapse=" + ")))
f1<-update(f0,.~.+map_q:period)
a<-anova(coxph(f0,data=obs),coxph(f1,data=obs),test="Chisq")
ints[["period"]]<-data.table(variable="period_observed_hours",p.value=a[2,"Pr(>|Chi|)"],n=nrow(obs),status="estimated")
it<-rbindlist(ints);it[,p_BH:=p.adjust(p.value,"BH")];fwrite(it,"02_Results/interactions.csv")

tfits<-lapply(dl,function(d){
 z<-survSplit(Surv(time_lm_days,event_lm)~.,data=as.data.frame(d[,c("time_lm_days","event_lm","map_q",primary),with=FALSE]),cut=6,episode="window")
 for(q in c("Q1","Q3","Q4"))z[[paste0(q,"_late")]]<-as.integer(z$map_q==q & z$window==2)
 f<-as.formula(paste("Surv(tstart,time_lm_days,event_lm) ~ map_q + Q1_late + Q3_late + Q4_late + strata(window) +",paste(primary,collapse=" + ")))
 coxph(f,data=z)
})
tt<-rbindlist(lapply(c(FALSE,TRUE),function(late){
 b<-vapply(tfits,function(f)coef(f)["map_qQ4"]+if(late)coef(f)["Q4_late"]else 0,numeric(1))
 v<-vapply(tfits,function(f){vv<-vcov(f);vv["map_qQ4","map_qQ4"]+if(late)vv["Q4_late","Q4_late"]+2*vv["map_qQ4","Q4_late"]else 0},numeric(1))
 z<-pool.scalar(b,v,n=tfits[[1]]$nevent,k=length(coef(tfits[[1]])))
 data.table(window=if(late)"ICU days 7-30"else"ICU days 1-7",HR=exp(z$qbar),lower=exp(z$qbar-qt(.975,z$df)*sqrt(z$t)),upper=exp(z$qbar+qt(.975,z$df)*sqrt(z$t)))
}))
fwrite(tt,"02_Results/time_specific_Q4.csv")
imp<-readRDS("restricted_cache/mice40_primary.rds")
cm<-as.data.table(as.table(imp$chainMean));setnames(cm,c("hour","iteration","chain","mean"));cm[,`:=`(iteration=as.integer(as.character(iteration)),chain=as.integer(gsub("[^0-9]","",as.character(chain))))]
fwrite(cm,"04_QC/chain_means.csv")
diags<-rbindlist(lapply(sp$hours,function(h){o<-imp$data[[h]][!is.na(imp$data[[h]])];a<-unlist(imp$imp[[h]]);data.table(hour=h,observed_n=length(o),missing_n=sum(is.na(imp$data[[h]])),observed_median=median(o),imputed_median=median(a),observed_q01=quantile(o,.01),observed_q99=quantile(o,.99),imputed_q01=quantile(a,.01),imputed_q99=quantile(a,.99),imputed_min=min(a),imputed_max=max(a))}))
fwrite(diags,"02_Results/imputation_distribution.csv")
conv<-tryCatch(as.data.table(mice::convergence(imp)),error=function(e)data.table(error=conditionMessage(e)))
fwrite(conv,"04_QC/MI_convergence.csv")
cat("Splines and diagnostics complete\n")
