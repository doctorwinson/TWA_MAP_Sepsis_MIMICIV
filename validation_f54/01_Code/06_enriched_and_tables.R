invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(data.table);library(mice);library(survival)})
sp<-readRDS("restricted_cache/model_spec.rds");asp<-readRDS("restricted_cache/analysis_specs.rds")
dl<-readRDS("restricted_cache/analysis_m40.rds");d<-dl[[1]]
imp<-readRDS("restricted_cache/mice40_enriched.rds")
edl<-lapply(1:40,function(i){z<-as.data.table(complete(imp,i));a<-copy(asp$raw);idx<-match(a$stay_id,z$stay_id);for(h in sp$hours)set(a,j=h,value=z[[h]][idx]);a[,mean_map:=rowMeans(.SD),.SDcols=sp$hours];a[,map_q:=relevel(cut(mean_map,c(-Inf,sp$cutpoints[2:4],Inf),labels=c("Q1","Q2","Q3","Q4"),right=TRUE),"Q2")];a})
specs<-list("Enriched MI primary"=asp$primary,"Enriched MI phenotype"=asp$pheno,"Enriched MI extended"=asp$extended,"Enriched MI SOFA components"=c(setdiff(asp$extended,"sofa"),asp$domains))
en<-rbindlist(lapply(names(specs),function(nm){
 fs<-lapply(edl,function(a)coxph(as.formula(paste("Surv(time_lm_days,event_lm) ~ map_q +",paste(specs[[nm]],collapse=" + "))),data=droplevels(a)))
 stopifnot(all(vapply(fs,function(f)all(is.finite(coef(f))),logical(1))))
 z<-as.data.table(summary(pool(fs),conf.int=TRUE));z[,`:=`(model=nm,HR=exp(estimate),lower=exp(conf.low),upper=exp(conf.high),n=fs[[1]]$n,events=fs[[1]]$nevent,m=40L)];z
}))
fwrite(en,"02_Results/enriched_MI_models.csv")

numeric_vars<-c("age","weight","charlson","sofa","mean_map","ttr65","ne_equiv_mean_0_24h","lactate",asp$domains)
binary_vars<-c("itvtn_24h_vent_tag","itvtn_24h_rrt_tag","drug_24h_sedative_tag","vaso_any","icd_hypertension","icd_stroke","icd_cad","cardio_proc_device","advanced_monitor","itvtn_24h_iabp_tag","itvtn_24h_cabg_tag","itvtn_24h_pci_tag","event_lm")
d[,female:=as.integer(gender=="F")];binary_vars<-c("female",binary_vars)
groups<-c("Overall","Q1","Q2","Q3","Q4")
labels<-c(age="Age, years",female="Female sex, n (%)",weight="Weight, kg",charlson="Charlson comorbidity index",sofa="SOFA score",mean_map="24-h TWA-MAP, mmHg",ttr65="Fraction of hours with MAP <65 mmHg",ne_equiv_mean_0_24h="NE-equivalent mean dose, micrograms/kg/min",lactate="First-day lactate, mmol/L",itvtn_24h_vent_tag="Mechanical ventilation, n (%)",itvtn_24h_rrt_tag="RRT, n (%)",drug_24h_sedative_tag="Sedative exposure, n (%)",vaso_any="Vasoactive infusion, n (%)",icd_hypertension="Hypertension code, n (%)",icd_stroke="Stroke code, n (%)",icd_cad="Coronary disease code, n (%)",cardio_proc_device="CABG/PCI/IABP, n (%)",advanced_monitor="PICCO/NICOM, n (%)",itvtn_24h_iabp_tag="IABP, n (%)",itvtn_24h_cabg_tag="CABG, n (%)",itvtn_24h_pci_tag="PCI, n (%)",event_lm="Death by ICU day 30, n (%)")
rows<-list();numeric_summary<-list()
for(v in c("age","female",setdiff(numeric_vars,"age"),setdiff(binary_vars,"female"))){
 lbl<-if(v %in% names(labels))labels[[v]]else gsub("_"," ",v)
 z<-data.table(variable=v,label=lbl)
 for(g in groups){
  a<-if(g=="Overall")d else d[map_q==g];x<-a[[v]]
  if(v %in% binary_vars){valid<-sum(!is.na(x));val<-sprintf("%d (%.1f%%)",sum(x==1,na.rm=TRUE),100*mean(x==1,na.rm=TRUE))}else{
   q<-as.numeric(quantile(x,c(.25,.5,.75),na.rm=TRUE));dig<-if(v=="age"||v %in% c("charlson","sofa",asp$domains))0 else if(v %in% c("ttr65","ne_equiv_mean_0_24h"))3 else 2
   fmt<-paste0("%.",dig,"f [%.",dig,"f, %.",dig,"f]");val<-sprintf(fmt,q[2],q[1],q[3]);valid<-sum(!is.na(x))
   numeric_summary[[paste(v,g)]]<-data.table(variable=v,group=g,n=valid,median=q[2],q25=q[1],q75=q[3])
  }
  z[,(g):=val];z[,(paste0(g,"_n")):=valid]
 }
 ref<-d[map_q=="Q2",get(v)]
 smds<-vapply(c("Q1","Q3","Q4"),function(g){a<-d[map_q==g,get(v)];a<-a[!is.na(a)];r<-ref[!is.na(ref)];den<-if(v %in% binary_vars)sqrt((mean(a)*(1-mean(a))+mean(r)*(1-mean(r)))/2)else sqrt((var(a)+var(r))/2);if(den==0)NA_real_ else abs(mean(a)-mean(r))/den},numeric(1))
 z[,max_abs_SMD_vs_Q2:=max(smds,na.rm=TRUE)];rows[[v]]<-z
}
fwrite(rbindlist(rows),"02_Results/table1_full.csv")
fwrite(rbindlist(numeric_summary),"02_Results/numeric_summaries.csv")
mainvars<-c("age","female","weight","charlson","sofa","mean_map","ttr65","ne_equiv_mean_0_24h","lactate","itvtn_24h_vent_tag","itvtn_24h_rrt_tag","drug_24h_sedative_tag","icd_hypertension","icd_stroke","event_lm")
fwrite(rbindlist(rows)[match(mainvars,variable)],"02_Results/table1_main.csv")
fwrite(data.table(variable=c(asp$primary,"lactate",asp$domains),missing=c(vapply(asp$raw[,asp$primary,with=FALSE],function(x)sum(is.na(x)),integer(1)),sum(is.na(asp$raw$lactate)),vapply(paste0(asp$domains,"_missing"),function(v)sum(asp$raw[[v]]),integer(1))),n=nrow(d)),"02_Results/covariate_missingness.csv")
desc<-data.table(metric=c("N","Deaths","complete24_N","median_missing_hours","q25_missing_hours","q75_missing_hours","observed_MAP_median","observed_MAP_q25","observed_MAP_q75","MAP_TTR_spearman","Q4_any_hypotension_percent","Q4_ttr25_percent"),
 value=c(nrow(d),sum(d$event_lm),sum(d$n_hour_obs==24),median(24-d$n_hour_obs),quantile(24-d$n_hour_obs,.25),quantile(24-d$n_hour_obs,.75),median(d$mean_observed),quantile(d$mean_observed,.25),quantile(d$mean_observed,.75),cor(d$mean_map,d$ttr65,method="spearman"),100*mean(d[map_q=="Q4",ttr65]>0),100*mean(d[map_q=="Q4",ttr65]>=.25)))
fwrite(desc,"02_Results/descriptive_metrics.csv")
cat("Enriched MI and editable table sources complete\n")
