invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(data.table);library(mice);library(survival);library(splines)})
setDTthreads(2)
logcon<-file("04_QC/models.log","wt");sink(logcon,split=TRUE)
sp<-readRDS("restricted_cache/model_spec.rds")
raw<-as.data.table(readRDS("restricted_cache/eligible_primary.rds"))
sv<-as.data.table(readRDS("restricted_cache/sepsis_source_verification.rds"))
if(!"both_pair_by24" %in% names(raw))raw<-merge(raw,sv[,.(stay_id,both_pair_by24)],by="stay_id",all.x=TRUE)
primary<-sp$primary;hours<-sp$hours
dom<-paste0("sofa_",c("respiration","coagulation","liver","cardiovascular","cns","renal"))
pheno<-c(primary,"icd_hypertension","icd_stroke","careunit_group","admission_group","cardio_proc_device")
extended<-c(pheno,"advanced_monitor","ne_log1p","service_group")
for(v in dom)raw[,(paste0(v,"_missing")):=as.integer(is.na(get(v)))]
raw[,complete_domains:=complete.cases(.SD),.SDcols=dom]
raw[,gender:=factor(gender)]
raw[,vaso_any:=as.integer(ne_equiv_mean_0_24h>0)]
raw[,lactate:=lab_fst_24h_lactate_first]
for(v in dom)set(raw,which(is.na(raw[[v]])),v,0)
stopifnot(all(rowSums(raw[,..dom])==raw$sofa))
raw[,sofa_noncv:=rowSums(.SD),.SDcols=setdiff(dom,"sofa_cardiovascular")]
raw[,sofa_noncns:=rowSums(.SD),.SDcols=setdiff(dom,"sofa_cns")]
keep<-unique(c("stay_id","subject_id","time_lm_days","event_lm",extended,dom,paste0(dom,"_missing"),"complete_domains","sofa_noncv","sofa_noncns","sapsii","vaso_any","both_pair_by24","period","n_hour_obs","mean_observed","lactate","hospital_expire_flag","ne_equiv_mean_0_24h","icd_cad","itvtn_24h_iabp_tag","itvtn_24h_cabg_tag","itvtn_24h_pci_tag",hours))
raw<-raw[,..keep]
categorize<-function(d) {
 d[,map_q:=relevel(cut(mean_map,c(-Inf,sp$cutpoints[2:4],Inf),labels=c("Q1","Q2","Q3","Q4"),right=TRUE),"Q2")]
 d
}
complete_list<-function(imp)lapply(1:40,function(i){
 z<-as.data.table(complete(imp,i))
 d<-copy(raw);idx<-match(d$stay_id,z$stay_id);stopifnot(!anyNA(idx))
 for(h in hours)set(d,j=h,value=z[[h]][idx])
 d[,mean_map:=rowMeans(.SD),.SDcols=hours]
 d[,ttr65:=rowMeans(.SD<65),.SDcols=hours]
 categorize(d)
})
dl<-complete_list(readRDS("restricted_cache/mice40_primary.rds"))
saveRDS(dl,"restricted_cache/analysis_m40.rds")
specs<-list("Crude"=character(),"Demographic"=primary[1:3],"Severity"=primary[1:5],"Primary"=primary,"Phenotype"=pheno,"Extended context"=extended,
 "Without added hospitalization codes"=c(primary,"careunit_group","admission_group","service_group","itvtn_24h_iabp_tag","advanced_monitor","ne_log1p"),
 "Plus vascular history"=c(primary,"icd_hypertension","icd_stroke"),"Plus ICU type"=c(primary,"careunit_group"),"Plus admission type"=c(primary,"admission_group"),"Plus service"=c(primary,"service_group"),"Plus procedures/devices"=c(primary,"cardio_proc_device","advanced_monitor"),"Plus NE dose"=c(primary,"ne_log1p"),
 "SOFA components"=c(setdiff(extended,"sofa"),dom),"Noncardiovascular SOFA"=c(setdiff(extended,"sofa"),"sofa_noncv"),"Non-CNS SOFA"=c(setdiff(extended,"sofa"),"sofa_noncns"),"Extended without MV"=setdiff(extended,"itvtn_24h_vent_tag"),"Extended without RRT"=setdiff(extended,"itvtn_24h_rrt_tag"),"Extended without sedation"=setdiff(extended,"drug_24h_sedative_tag"),
 "SAPS II replacement"=c(setdiff(primary,"sofa"),"sapsii"),"Lactate and NE dose"=c(primary,"lactate","ne_log1p"))
fit_storage<-list();warns<-list()
run<-function(label,covariates,dds=dl,filter=NULL,exposure="map_q",save=FALSE){
 cat(label,"\n")
 ww<-character()
 fits<-tryCatch(withCallingHandlers(lapply(dds,function(d){
  if(!is.null(filter))d<-d[filter(d)]
  needed<-unique(c("time_lm_days","event_lm",exposure,covariates))
  d<-droplevels(as.data.frame(d[complete.cases(d[,..needed])]))
  cv<-covariates[vapply(d[,covariates,drop=FALSE],function(x)length(unique(x))>1,logical(1))]
  if(nrow(d)<100 || sum(d$event_lm)<25)stop("Sparse sample")
  f<-as.formula(paste("Surv(time_lm_days,event_lm) ~",paste(c(exposure,cv),collapse=" + ")))
  coxph(f,data=d,ties="efron",x=save,model=save)
 }),warning=function(w){ww<<-c(ww,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
 if(inherits(fits,"error"))return(data.table(model=label,term="map_qQ4",status=conditionMessage(fits)))
 if(any(vapply(fits,function(f)any(!is.finite(coef(f)))||any(!is.finite(vcov(f))),logical(1)))||any(grepl("infinite|converge",ww)))return(data.table(model=label,term="map_qQ4",n=fits[[1]]$n,events=fits[[1]]$nevent,m=length(fits),status=paste("unstable",paste(unique(ww),collapse=";"))))
 if(length(ww))warns[[label]]<<-unique(ww)
 if(save)fit_storage[[label]]<<-fits
 if(length(fits)>1) {
  out<-as.data.table(summary(pool(fits),conf.int=TRUE))
 }else{
  b<-coef(fits[[1]]);se<-sqrt(diag(vcov(fits[[1]])))
  out<-data.table(term=names(b),estimate=b,std.error=se,conf.low=b-1.96*se,conf.high=b+1.96*se,p.value=2*pnorm(-abs(b/se)))
 }
 out[,`:=`(model=label,HR=exp(estimate),lower=exp(conf.low),upper=exp(conf.high),n=fits[[1]]$n,events=fits[[1]]$nevent,m=length(fits),status="estimated")]
 out
}
models<-rbindlist(lapply(names(specs),function(nm)run(nm,specs[[nm]],save=nm %in% c("Primary","Phenotype","Extended context"))),fill=TRUE)
fwrite(models,"02_Results/main_models.csv")
q4<-models[term=="map_qQ4"]
q4[,descriptive_logHR_change_pct:=100*(q4[model=="Primary",estimate]-estimate)/q4[model=="Primary",estimate]]
fwrite(q4,"02_Results/phenotype_Q4.csv")
sens<-rbindlist(list(
 run("Observed >=18 hours",primary,list(categorize(copy(raw)[,mean_map:=mean_observed]))),
 run("Observed h02-h23",primary,list(categorize(copy(raw)[,mean_map:=rowMeans(.SD,na.rm=TRUE),.SDcols=hours[3:24]]))),
 run("Complete 24 hours",primary,list(categorize(copy(raw)[n_hour_obs==24,mean_map:=mean_observed][n_hour_obs==24]))),
 run("Complete 24 hours extended",extended,list(categorize(copy(raw)[n_hour_obs==24,mean_map:=mean_observed][n_hour_obs==24]))),
 run("At least 20 observed hours",primary,filter=function(d)d$n_hour_obs>=20),
 run("Antibiotic/culture pair completed by 24h",primary,filter=function(d)d$both_pair_by24),
 run("Pair completed by 24h extended",extended,filter=function(d)d$both_pair_by24),
 run("Complete raw SOFA components",c(setdiff(extended,"sofa"),dom),filter=function(d)d$complete_domains)
 ),fill=TRUE)
fwrite(sens,"02_Results/sensitivity_models.csv")
filters<-list("Medical ICU"=function(d)d$careunit_group=="Medical ICU","Cardiac ICU"=function(d)d$careunit_group=="Cardiac ICU","Surgical/Trauma ICU"=function(d)d$careunit_group=="Surgical/Trauma ICU","Neuro ICU/Stepdown"=function(d)d$careunit_group=="Neuro ICU/Stepdown","No MV"=function(d)d$itvtn_24h_vent_tag==0,"MV"=function(d)d$itvtn_24h_vent_tag==1,"No RRT"=function(d)d$itvtn_24h_rrt_tag==0,"RRT"=function(d)d$itvtn_24h_rrt_tag==1,"No vasoactive infusion"=function(d)d$vaso_any==0,"Vasoactive infusion"=function(d)d$vaso_any==1,"Exclude neuro ICU"=function(d)d$careunit_group!="Neuro ICU/Stepdown","Exclude cardiac procedures/devices"=function(d)d$cardio_proc_device==0)
subgroups<-rbindlist(lapply(names(filters),function(nm)run(nm,primary,filter=filters[[nm]])),fill=TRUE)
fwrite(subgroups,"02_Results/subgroup_models.csv")
tc<-list()
for(p in levels(raw$period))for(adj in c("Primary","Extended context")) {
 obs<-categorize(copy(raw)[period==p,mean_map:=mean_observed][period==p])
 tc[[paste(p,adj)]]<-run(paste(p,adj),specs[[adj]],list(obs))
}
fwrite(rbindlist(tc,fill=TRUE),"02_Results/temporal_models.csv")
fwrite(dl[[1]][,.(n=.N,events=sum(event_lm),mortality_pct=100*mean(event_lm)),by=map_q][order(map_q)],"02_Results/mortality_by_quartile.csv")
fwrite(raw[,.(n=.N,events=sum(event_lm)),by=period],"02_Results/temporal_counts.csv")
catvars<-c("careunit_group","service_group","admission_group","cardio_proc_device","advanced_monitor","itvtn_24h_vent_tag","itvtn_24h_rrt_tag","drug_24h_sedative_tag","vaso_any","icd_hypertension","icd_stroke","icd_cad","itvtn_24h_iabp_tag","itvtn_24h_cabg_tag")
desc<-rbindlist(lapply(catvars,function(v){x<-dl[[1]][,.N,by=c("map_q",v)];setnames(x,c(v,"N"),c("level","n"));x[,level:=as.character(level)];x[,denominator:=sum(n),by=map_q];x[,`:=`(variable=v,percent=100*n/denominator)];x}))
fwrite(desc,"02_Results/phenotype_distribution.csv")
saveRDS(list(primary=primary,pheno=pheno,extended=extended,domains=dom,raw=raw),"restricted_cache/analysis_specs.rds")
saveRDS(fit_storage,"restricted_cache/core_fits.rds")
if(length(warns))capture.output(warns,file="04_QC/model_warnings.txt")
capture.output(sessionInfo(),file="04_QC/model_session.txt")
sink();close(logcon)
