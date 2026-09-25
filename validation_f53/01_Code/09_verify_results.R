invisible(Sys.setlocale("LC_CTYPE","English_United States.utf8"))
suppressPackageStartupMessages({library(data.table);library(mice)})
checks<-list()
check<-function(name,ok,detail=""){checks[[name]]<<-data.table(check=name,passed=isTRUE(ok),detail=detail)}
files<-list.files("01_Code",pattern="\\.R$",full.names=TRUE)
for(f in files)check(paste("R parse",basename(f)),!inherits(try(parse(f),silent=TRUE),"try-error"))
mi<-readRDS("restricted_cache/mice40_full.rds")$imp
check("m=40 and ten iterations",mi$m==40 && mi$iteration==10)
check("MI population 7766",nrow(mi$data)==7766)
hc<-sprintf("h%02d",0:23)
for(i in 1:40){d<-complete(mi,i);check(paste("Completed hour matrices",i),!anyNA(d[,hc]))}
r<-fread("02_Results/P00_F52_reproduction.csv")
p<-r[model=="Primary all eligible" & term=="map_qQ4"]
q<-r[model=="Phenotype all eligible" & term=="map_qQ4"]
check("F52 primary reproduced",round(p$HR,2)==2.18 && round(p$lower,2)==1.82 && round(p$upper,2)==2.62)
check("F52 phenotype reproduced",round(q$HR,2)==1.55 && round(q$lower,2)==1.28 && round(q$upper,2)==1.87)
cuts<-fread("02_Results/rebuilt_cutpoints.csv")
check("F52 cutpoints reproduced",all(round(cuts$value[2:4],2)==c(71.40,75.72,80.99)))
b<-fread("02_Results/P08_full_cohort_Q4_changes.csv")
check("Full cohort block comparisons use identical n",all(b$n==7765) && all(b$events==1069))
ee<-fread("02_Results/E03_harmonized_hospital_models.csv")[model=="eICU shared primary" & term=="map_qQ4"]
check("eICU denominator and endpoint",ee$n==2765 && ee$events==714 && ee$hospitals==153 && ee$endpoint=="In-hospital mortality")
check("eICU OR interval finite",is.finite(ee$OR)&&ee$lower<ee$OR&&ee$OR<ee$upper)
int<-fread("02_Results/V06_pooled_interactions.csv")
check("Six exploratory interactions retained",nrow(int)==6 && !anyNA(int$p.value) && all(int$p_BH>=int$p.value-1e-12))
raw<-readRDS("restricted_cache/mimic_enriched_preimputation.rds")
check("Timing audit 764 late sepsis",sum(raw$onset_hours>24)==764)
for(f in c("mice40_enriched_full.rds","mice40_enriched_complete_domains.rds")) {
 m<-readRDS(file.path("restricted_cache",f))$imp
 check(paste("Enriched MI",f),m$m==40L&&m$iteration==10L)
 if(!is.null(m$loggedEvents)) {
  z<-as.data.table(m$loggedEvents)[,.N,by=.(meth,out)]
  fwrite(z,file.path("04_QC",paste0(f,"_events_summary.csv")))
 }
}
for(f in list.files("02_Results",pattern="\\.csv$",full.names=TRUE,recursive=TRUE)) {
 cols<-names(fread(f,nrows=0))
 check(paste("Aggregate-only columns",basename(f)),!any(cols %in% c("stay_id","subject_id","hadm_id","uniquepid","patientunitstayid","hospitalid")))
}
out<-rbindlist(checks);fwrite(out,"04_QC/numerical_verification.csv")
cat(sum(out$passed),"/",nrow(out),"checks passed\n")
if(any(!out$passed))stop("Numerical verification failed")
